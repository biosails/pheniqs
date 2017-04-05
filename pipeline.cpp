/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pipeline.h"

/*  Pivot */

Pivot::Pivot(Pipeline& pipeline) :
    index(pipeline.pivots.size()),
    decoder(pipeline.environment.decoder),
    multiplex_barcode(pipeline.environment.total_multiplex_barcode_segments),
    molecular_barcode(pipeline.environment.total_molecular_barcode_segments),
    accumulator(pipeline.environment.total_input_segments),
    pipeline(pipeline),
    disable_quality_control(pipeline.environment.disable_quality_control),
    long_read(pipeline.environment.long_read),
    leading_segment(NULL) {

    // create input segments
    input.reserve(pipeline.environment.total_input_segments);
    for(size_t i = 0; i < pipeline.environment.total_input_segments; i++) {
        input.emplace_back(i, i + 1, pipeline.environment.total_input_segments, pipeline.environment.platform);
    }

    // create output segments
    output.reserve(pipeline.environment.total_output_segments);
    for(size_t i = 0; i < pipeline.environment.total_output_segments; i++) {
        output.emplace_back(i, i + 1, pipeline.environment.total_output_segments, pipeline.environment.platform);
    }

    // set first output segment READ1 flag ON
    if(pipeline.environment.total_output_segments > 0) {
        output[0].flag |= uint16_t(HtsFlag::READ1);
    }

    // set last output segment READ2 flag ON
    if(pipeline.environment.total_output_segments > 1) {
        output[pipeline.environment.total_output_segments - 1].flag |= uint16_t(HtsFlag::READ2);
    }

    leading_segment = &(input[pipeline.environment.leading_segment_index]);

    /*
        in long read accumulation more channel statistics
        are collected on the pivots and only at the end summed up
        on the channel. otherwise they are collected directly on the channel
    */
    if(long_read) {
        for(const auto specification : pipeline.environment.channel_specifications) {
            channel_by_barcode.emplace(make_pair(string(specification->multiplex_barcode), ChannelAccumulator(*specification)));
        }
    }
    clear();
};
void Pivot::start() {
    pivot_thread = thread(&Pivot::run, this);
};
inline void Pivot::clear() {
    // determined = false;
    // filtered = false;
    // multiplex_probability = 0;
    // conditioned_multiplex_probability = 0;
    // multiplex_distance = 0;
    // decoded_multiplex_channel = NULL;
    multiplex_barcode.clear();
    molecular_barcode.clear();
    for(auto& segment : input) segment.clear();
    for(auto& segment : output) segment.clear();
};
void Pivot::run() {
    switch (pipeline.environment.action) {
    case ProgramAction::DEMULTIPLEX: {
        switch(decoder) {
        case(Decoder::PAMLD): {
            while(pipeline.pull(*this)) {
                validate();
                transform();
                decode_with_pamld();
                encode_auxiliary();
                push();
                increment();
                clear();
            }
            break;
        };
        case(Decoder::MDD): {
            while(pipeline.pull(*this)) {
                validate();
                transform();
                decode_with_mdd();
                encode_auxiliary();
                push();
                increment();
                clear();
            }
            break;
        };
        case(Decoder::BENCHMARK): {
            while(pipeline.pull(*this)) {
                validate();
                transform();
                decode_with_mdd();
                encode_mmd_benchmark_auxiliary();
                decode_with_pamld();
                encode_pamld_benchmark_auxiliary();
                encode_auxiliary();
                push();
                increment();
                clear();
            }
            break;
        };
        default: break;
        }
        break;
    };
    case ProgramAction::QUALITY: {
        while(pipeline.pull(*this)) {
            validate();
            transform();
            decode_tagged_channel();
            copy_auxiliary();
            push();
            increment();
            clear();
        }
        break;
    };
    default: break;
    }
};
void Pivot::rotate() {
    while(pipeline.pull(*this)) {
        transform();
        copy_auxiliary();
        if(leading_segment->auxiliary.RG.l > 0) {
            string rg(leading_segment->auxiliary.RG.s, leading_segment->auxiliary.RG.l);
            auto record = pipeline.channel_by_read_group_id.find(rg);
            if (record != pipeline.channel_by_read_group_id.end()) {
                record->second->push(*this);
            }
        }
        increment();
    }
};
inline void Pivot::validate() {
    if(input.size() > 1) {
        // validate that all segments in the pivot have the same identifier
        const kstring_t& baseline = input[0].name;
        for(size_t i = 1; i < input.size(); i++) {
            Segment& segment = input[i];
            if((baseline.l != segment.name.l) || strncmp(baseline.s, segment.name.s, baseline.l)) {
                throw SequenceError("read segments out of sync " + string(segment.name.s, segment.name.l) + " and " + string(baseline.s, baseline.l));
            }
        }
    }
};
inline void Pivot::decode_with_pamld() {
    multiplex_distance = pipeline.environment.multiplex_barcode_width();
    decoded_multiplex_channel = pipeline.undetermined;
    conditioned_multiplex_probability = 0;
    multiplex_probability = 0;
    determined = false;

    /*  Compute P(observed|barcode) for each barcode
        Keep track of the channel that yield the maximal prior adjusted probability.
        If r is the observed sequence and b is the barcode sequence
        P(r|b) is the probability that r was observed given b was sequenced.
        Accumulate all prior adjusted probabilities P(b) * P(r|b), in sigma
        using the Kahan summation algorithm to minimize floating point drift
        see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    */
    double adjusted = 0;
    double compensation = 0;
    double sigma = 0;
    double y = 0;
    double t = 0;
    double c = 0;
    double p = 0;
    size_t d = 0;
    for (auto channel : pipeline.channels) {
        channel->multiplex_barcode.accurate_decoding_probability(multiplex_barcode, c, d);
        p = c * channel->concentration;
        y = p - compensation;
        t = sigma + y;
        compensation = (t - sigma) - y;
        sigma = t;
        if (p > adjusted) {
            decoded_multiplex_channel = channel;
            conditioned_multiplex_probability = c;
            multiplex_distance = d;
            adjusted = p;
        }
    }

    /*  Compute P(barcode|observed)
        P(b|r), the probability that b was sequenced given r was observed
        P(b|r) = P(r|b) * P(b) / ( P(noise) * P(r|noise) + sigma )
        where sigma = sum of P(r|b) * P(b) over b
    */
    multiplex_probability = adjusted / (sigma + pipeline.environment.adjusted_noise_probability);

    /* Check for decoding failure and assign to the undetermined channel if decoding failed */
    if (conditioned_multiplex_probability > pipeline.conditioned_probability_threshold &&
        multiplex_probability > pipeline.environment.confidence) {
        determined = true;

    } else {
        decoded_multiplex_channel = pipeline.undetermined;
        conditioned_multiplex_probability = 0;
        multiplex_distance = 0;
        multiplex_probability = 1;
        // multiplex_probability = 1.0 - sigma / (sigma + pipeline.environment.adjusted_noise_probability);
    }
};
inline void Pivot::decode_with_mdd() {
    multiplex_distance = pipeline.environment.multiplex_barcode_width();
    decoded_multiplex_channel = pipeline.undetermined;
    determined = false;

    /* First try a perfect match to the full barcode sequence */
    const string code(multiplex_barcode);
    auto record = pipeline.channel_by_barcode.find(code);
    if (record != pipeline.channel_by_barcode.end()) {
        decoded_multiplex_channel = record->second;
        multiplex_distance = 0;
        determined = true;

    } else {
        /* If no exact match was found try error correction */
        for (auto channel : pipeline.channels) {
            if (channel->multiplex_barcode.corrected_match(multiplex_barcode, multiplex_distance)) {
                decoded_multiplex_channel = channel;
                determined = true;
                break;
            }
        }
        if(!determined) {
            multiplex_distance = 0;
        }
    }
};
inline void Pivot::transform() {
    // populate the output
    for (auto& transform : pipeline.environment.template_transforms ) {
        const Segment& from = input[transform.token.input_segment_index];
        Segment& to = output[transform.output_segment_index];
        to.sequence.append(from.sequence, transform);
    }

    // populate the multiplex barcode
    for (auto& transform : pipeline.environment.multiplex_barcode_transforms) {
        const Segment& from = input[transform.token.input_segment_index];
        multiplex_barcode.append(transform.output_segment_index, from.sequence, transform);
    }

    // populate the molecular barcode
    for (auto& transform : pipeline.environment.molecular_barcode_transforms) {
        const Segment& from = input[transform.token.input_segment_index];
        molecular_barcode.append(transform.output_segment_index, from.sequence, transform);
    }

    // assign the pivot qc_fail flag from the leader
    filtered = leading_segment->get_qcfail();

    for (auto& segment : output) {
        kputsn(leading_segment->name.s, leading_segment->name.l, &segment.name);
        segment.set_qcfail(leading_segment->get_qcfail());
        segment.auxiliary.XI = leading_segment->auxiliary.XI;
    }
};
inline void Pivot::decode_tagged_channel() {
    if(leading_segment->auxiliary.RG.l > 0) {
        string rg(leading_segment->auxiliary.RG.s, leading_segment->auxiliary.RG.l);
        auto record = pipeline.channel_by_read_group_id.find(rg);
        if (record != pipeline.channel_by_read_group_id.end()) {
            decoded_multiplex_channel = record->second;
        }
    }
};
inline void Pivot::encode_auxiliary () {
    if (decoded_multiplex_channel != NULL) {
        for(auto& segment: output) {

            if(decoder == Decoder::PAMLD) 
                segment.auxiliary.DQ = 1.0 - multiplex_probability;

            // sequence auxiliary tags
            segment.sequence.expected_error(segment.auxiliary.EE);

            // pivot auxiliary tags
            segment.auxiliary.set_multiplex_barcode(multiplex_barcode);
            segment.auxiliary.set_molecular_barcode(molecular_barcode);

            // channel auxiliary tags
            kputsn(decoded_multiplex_channel->rg.ID.s, decoded_multiplex_channel->rg.ID.l, &segment.auxiliary.RG);
         /* kputsn(decoded_multiplex_channel->rg.LB.s, decoded_multiplex_channel->rg.LB.l, &segment.auxiliary.LB);
            kputsn(decoded_multiplex_channel->rg.PU.s, decoded_multiplex_channel->rg.PU.l, &segment.auxiliary.PU);
            kputsn(decoded_multiplex_channel->rg.FS.s, decoded_multiplex_channel->rg.FS.l, &segment.auxiliary.FS);
            kputsn(decoded_multiplex_channel->rg.PG.s, decoded_multiplex_channel->rg.PG.l, &segment.auxiliary.PG);
            kputsn(decoded_multiplex_channel->rg.CO.s, decoded_multiplex_channel->rg.CO.l, &segment.auxiliary.CO); */
        }
    }
};
inline void Pivot::copy_auxiliary () {
    for(auto& segment: output) {
        segment.auxiliary.imitate(leading_segment->auxiliary);
    }
};
inline void Pivot::push() {
    if (decoded_multiplex_channel != NULL) {
        decoded_multiplex_channel->push(*this);
    }
};
inline void Pivot::increment() {
    if(!disable_quality_control) {
        // input quality tracking
        accumulator.increment(*this);

        // long read tracks channel quality on the pivot
        if(long_read) {
            const auto& record = channel_by_barcode.find(decoded_multiplex_channel->key);
            if (record != channel_by_barcode.end()) {
                record->second.increment(*this);
            }
        }
    }
};
inline void Pivot::encode_mmd_benchmark_auxiliary () {
    if(decoded_multiplex_channel != NULL) {
        for(auto& segment: output) {
            segment.auxiliary.set_mdd_barcode(decoded_multiplex_channel->multiplex_barcode);
            segment.auxiliary.YD = multiplex_distance;
        }
    }
};
inline void Pivot::encode_pamld_benchmark_auxiliary () {
    if(decoded_multiplex_channel != NULL) {
        for(auto& segment: output) {
            segment.auxiliary.set_pamld_barcode(decoded_multiplex_channel->multiplex_barcode);
            segment.auxiliary.XD = multiplex_distance;
            segment.auxiliary.XP = conditioned_multiplex_probability;
        }
    }
};

/*  Channel */

Channel::Channel(const Pipeline& pipeline, const ChannelSpecification& specification) :
    index(specification.index),
    key(specification.multiplex_barcode),
    concentration(specification.concentration),
    multiplex_barcode(specification.multiplex_barcode),
    disable_quality_control(specification.disable_quality_control),
    long_read(specification.long_read),
    include_filtered(specification.include_filtered),
    undetermined(specification.undetermined),
    rg(specification.rg),
    channel_accumulator(specification) {

    map< size_t, Feed* > feed_by_index;
    output_feeds.reserve(specification.output_urls.size());
    for(const auto& url : specification.output_urls) {
        Feed* feed = pipeline.resolve_feed(url, IoDirection::OUT);
        output_feeds.push_back(feed);
        if (feed_by_index.find(feed->index()) == feed_by_index.end()) {
            feed_by_index[feed->index()] = feed;
        }
    }
    unique_output_feeds.reserve(feed_by_index.size());
    for(const auto& record : feed_by_index) {
        unique_output_feeds.push_back(record.second);
    }
};
Channel::~Channel() {
};
inline void Channel::increment(Pivot& pivot) {
    lock_guard<mutex> channel_lock(channel_mutex);
    channel_accumulator.increment(pivot);
};
void Channel::push(Pivot& pivot) {
    // acquire a push lock for all feeds in a fixed order
    vector< unique_lock< mutex > > feed_locks;
    feed_locks.reserve(unique_output_feeds.size());
    for(const auto feed : unique_output_feeds) {
        feed_locks.push_back(feed->acquire_push_lock());
    }

    // push the segments to the output feeds
    for(size_t i = 0; i < output_feeds.size(); i++) {
        output_feeds[i]->push(pivot.output[i]);
    }

    // release the locks on the feeds in reverse order
    for (auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }

    // for short reads the channels track their own quality
    if(!disable_quality_control && !long_read) {
        increment(pivot);
    }
};
void Channel::encode(Document& document, Value& value, const bool disable_quality_control) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    Value v;

    v.SetUint64(index);
    value.AddMember("index", v, allocator);

    // if (!multiplex_barcode.empty()) {
    //     Value barcodes;
    //     barcodes.SetArray();
    //     for (auto& sequence : multiplex_barcode.fragments) {
    //         Value barcode;
    //         barcode.SetObject();
    //         string code = sequence.to_iupac_ambiguity();
    //         v.SetString(code.c_str(), code.size(), allocator);
    //         barcode.AddMember("barcode sequence", v, allocator);
    //         barcodes.PushBack(barcode, allocator);
    //     }
    //     value.AddMember("barcodes", barcodes, allocator);
    // }
    channel_accumulator.encode(document, value, disable_quality_control);
};
void Channel::finalize() {
    channel_accumulator.finalize();
};

/*  Pipeline */

Pipeline::Pipeline(Environment& environment) :
    environment(environment),
    conditioned_probability_threshold(environment.random_word_probability),
    disable_quality_control(environment.disable_quality_control),
    long_read(environment.long_read),
    end_of_input(false),
    thread_pool({NULL, 0}),
    undetermined(NULL),
    input_accumulator(environment.total_input_segments) {

    thread_pool.pool = hts_tpool_init(environment.threads);
    if (!thread_pool.pool) {
        throw InternalError("error creating thread pool");
    }

    environment.validate_urls();
    load_input_feeds();
    load_output_feeds();
};
Pipeline::~Pipeline() {
    hts_tpool_destroy(thread_pool.pool);
    for(auto feed : unique_input_feeds) {
        delete feed;
    }
    for(auto feed : unique_output_feeds) {
        delete feed;
    }
    input_feeds.clear();
    unique_input_feeds.clear();
    unique_output_feeds.clear();

    for(auto pivot : pivots) {
        delete pivot;
    }
    for(auto channel : channels) {
        delete channel;
    }
    delete undetermined;
    undetermined = NULL;
    channels.clear();
    channel_by_barcode.clear();
    channel_by_read_group_id.clear();
};
void Pipeline::initialize_channels() {
    channels.reserve(environment.channel_specifications.size());
    for(const auto specification : environment.channel_specifications) {
        Channel* channel = new Channel(*this, *specification);
        if(!channel->undetermined) {
            channels.push_back(channel);

        } else {
            // the undetermined channel is not present in the channels array
            // and is referenced separately by the undetermined pointer
            undetermined = channel;
        }

        // channel by read group id lookup
        string id(channel->rg.ID.s, channel->rg.ID.l);
        if (channel_by_read_group_id.find(id) == channel_by_read_group_id.end()) {
            channel_by_read_group_id.emplace(make_pair(id, channel));
        }

        if (channel_by_barcode.find(channel->key) == channel_by_barcode.end()) {
            channel_by_barcode[channel->key] = channel;
        }
    }
};
void Pipeline::initialize_pivots() {
    pivots.reserve(environment.transforms);
    for(size_t i = 0; i < environment.transforms; i++) {
        Pivot* pivot = new Pivot(*this);
        pivots.push_back(pivot);
    }
};
void Pipeline::execute() {
    // initialize differently if quality
    // can already inspect the input
    initialize_channels();
    initialize_pivots();

    start();
    for (auto pivot : pivots) {
        pivot->pivot_thread.join();
    }
    stop();

    // consolidate accumulators from all pivots
    finalize();

    // write report to stdout
    encode_report(cout);
};
bool Pipeline::pull(Pivot& pivot) {
    vector< unique_lock< mutex > > feed_locks;
    feed_locks.reserve(unique_input_feeds.size());

    // acquire a pull lock for all feeds in a fixed order
    for(const auto feed : unique_input_feeds) {
        feed_locks.push_back(feed->acquire_pull_lock());
    }

    // pull into pivot input segments from input feeds
    for(size_t i = 0; i < pivot.input.size(); i++) {
        if(!input_feeds[i]->pull(pivot.input[i])) {
            end_of_input = true;
        }
    }

    // release the locks on the feeds in reverse order
    for (auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }
    return !end_of_input;
};
void Pipeline::encode_report(ostream& o) const {
    switch (environment.action) {
        case ProgramAction::DEMULTIPLEX: {
            encode_demultiplex_report(o);
            break;
        };
        case ProgramAction::QUALITY: {
            encode_quality_report(o);
            break;
        };
        default:
            break;
    }
};
void Pipeline::encode_input_report(Document& document, Value& value) const {
    value.SetObject();
    input_accumulator.encode(document, value, disable_quality_control);
};
void Pipeline::encode_output_report(Document& document, Value& value) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    value.SetObject();

    uint64_t total_yield = 0;
    uint64_t total_perfect_yield = 0;
    uint64_t total_filtered_yield = 0;

    for (auto channel : channels) {
        total_yield += channel->channel_accumulator.yield;
        total_perfect_yield += channel->channel_accumulator.perfect_yield;
        total_filtered_yield += channel->channel_accumulator.filtered_yield;
    }

    Value v;
    v.SetUint64(total_yield);
    value.AddMember("yield", v, allocator);

    v.SetUint64(total_perfect_yield);
    value.AddMember("perfect yield", v, allocator);

    v.SetUint64(total_filtered_yield);
    value.AddMember("filtered yield", v, allocator);

    v.SetDouble(total_yield > 0 ? double(total_perfect_yield) / double(total_yield) : 0.0);
    value.AddMember("perfect yield fraction", v, allocator);

    v.SetDouble(total_yield > 0 ? double(total_filtered_yield) / double(total_yield) : 0.0);
    value.AddMember("filtered yield fraction", v, allocator);

    Value channel_reports;
    channel_reports.SetArray();
    for (size_t i = 0; i < channels.size(); i++) {
        Channel* channel = channels[i];

        Value channel_report;
        channel_report.SetObject();
        channel->encode(document, channel_report, disable_quality_control);
        channel_reports.PushBack(channel_report, allocator);

        // Value segmentReports;
        // segmentReports.SetArray();
        // for (uint_fast32_t i = 0; i < channel->output_feeds.size(); i++) {
        //     Value segmentReport;
        //     segmentReport.SetObject();
        //     v.SetString((char*)channel->output_feeds[i]->url.c_str(), allocator);
        //     segmentReport.AddMember("path", v, allocator);
        //     if (!disable_quality_control) {
        //         channel->distribution.feed_statistics[i].encodeReport(document, &segmentReport);
        //     }
        //     segmentReports.PushBack(segmentReport, allocator);
        // }
    }
    value.AddMember("library quality reports", channel_reports, allocator);
};
void Pipeline::encode_demultiplex_report(ostream& o) const {
    Document document;
    document.SetObject();
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;

    v.SetUint64(input_accumulator.count);
    document.AddMember("read count", v, allocator);

    v.SetUint64(input_accumulator.filtered_count);
    document.AddMember("filtered read count", v, allocator);

    v.SetDouble(input_accumulator.count > 0 ? double(input_accumulator.filtered_count) / double(input_accumulator.count) : 0.0);
    document.AddMember("filtered read fraction", v, allocator);

    v.SetUint64(input_accumulator.determined_count);
    document.AddMember("determined read count", v, allocator);

    v.SetUint64(input_accumulator.determined_filtered_count);
    document.AddMember("determined filtered read count", v, allocator);

    v.SetDouble(input_accumulator.count > 0 ? double(input_accumulator.determined_count) / double(input_accumulator.count) : 0.0);
    document.AddMember("determined read fraction", v, allocator);

    Value output_report;
    encode_output_report(document, output_report);
    document.AddMember("demultiplex output report", output_report, allocator);

    Value input_report;
    encode_input_report(document, input_report);
    document.AddMember("demultiplex input report", input_report, allocator);

    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    document.Accept(writer);
    o << buffer.GetString();
};
void Pipeline::encode_quality_report(ostream& o) const {
    Document document;
    document.SetObject();
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;

    v.SetUint64(input_accumulator.count);
    document.AddMember("read count", v, allocator);

    v.SetUint64(input_accumulator.filtered_count);
    document.AddMember("filtered read count", v, allocator);

    v.SetDouble(input_accumulator.count > 0 ? double(input_accumulator.filtered_count) / double(input_accumulator.count) : 0.0);
    document.AddMember("filtered read fraction", v, allocator);

    Value output_report;
    encode_output_report(document, output_report);
    document.AddMember("demultiplex output report", output_report, allocator);

    Value input_report;
    encode_input_report(document, input_report);
    document.AddMember("demultiplex input report", input_report, allocator);

    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    document.Accept(writer);
    o << buffer.GetString();
};
void Pipeline::finalize() {
    // collect statistics from all parallel pivots
    for(const auto pivot : pivots) {
        input_accumulator += pivot->accumulator;

        if(long_read) {
            for(const auto& record : pivot->channel_by_barcode) {
                const auto& channel = channel_by_barcode.find(record.first);
                if (channel != channel_by_barcode.end()) {
                    channel->second->channel_accumulator += record.second;
                }
            }
        }
    }
    input_accumulator.finalize();

    for(auto channel : channels) {
        channel->finalize();
    }
};
Feed* Pipeline::resolve_feed(const URL& url, const IoDirection& direction) const {
    Feed* feed = NULL;
    switch(direction) {
        case IoDirection::IN: {
            const auto& record = input_feed_by_url.find(url);
            if(record != input_feed_by_url.end()) {
                feed = record->second;
            }
            break;
        };
        case IoDirection::OUT: {
            const auto& record = output_feed_by_url.find(url);
            if(record != output_feed_by_url.end()) {
                feed = record->second;
            }
            break;
        };
    }
    return feed;
};
void Pipeline::load_input_feeds() {
    map< size_t, Feed* > feed_by_index;

    input_feed_by_url.reserve(environment.input_feed_specification_by_url.size());
    for(const auto& record : environment.input_feed_specification_by_url) {
        FeedSpecification* specification = record.second;
        Feed* feed = NULL;

        switch(specification->url.kind()) {
            case FormatKind::FASTQ: {
                feed = new FastqFeed(specification);
                break;
            };
            case FormatKind::HTS: {
                feed = new HtsFeed(specification);
                break;
            };
            default:
                throw InternalError("unknown input format " + specification->url);
                break;
        }
        input_feed_by_url[specification->url] = feed;
        feed_by_index[feed->index()] = feed;
    }

    input_feeds.reserve(environment.input_urls.size());
    for(const auto& url : environment.input_urls) {
        Feed* feed = resolve_feed(url, IoDirection::IN);
        input_feeds.push_back(feed);
    }

    unique_input_feeds.reserve(feed_by_index.size());
    for(const auto& record : feed_by_index) {
        unique_input_feeds.push_back(record.second);
    }
    for(auto feed : unique_input_feeds) {
        feed->set_thread_pool(&thread_pool);
    }
};
void Pipeline::load_output_feeds() {
    map< size_t, Feed* > feed_by_index;

    output_feed_by_url.reserve(environment.output_feed_specification_by_url.size());
    for(const auto& record : environment.output_feed_specification_by_url) {
        FeedSpecification* specification = record.second;
        Feed* feed = NULL;

        switch(specification->url.kind()) {
            case FormatKind::FASTQ: {
                feed = new FastqFeed(specification);
                break;
            };
            case FormatKind::HTS: {
                feed = new HtsFeed(specification);
                break;
            };
            default:
                throw InternalError("unknown output format " + specification->url);
                break;
        }
        output_feed_by_url[specification->url] = feed;
        feed_by_index[feed->index()] = feed;
    }

    unique_output_feeds.reserve(feed_by_index.size());
    for(const auto& record : feed_by_index) {
        unique_output_feeds.push_back(record.second);
    }
    for(auto feed : unique_output_feeds) {
        feed->set_thread_pool(&thread_pool);
    }
};
void Pipeline::start() {
    for(auto feed : unique_input_feeds) {
        feed->open();
    }
    for(auto feed : unique_output_feeds) {
        feed->open();
    }
    for(auto feed : unique_input_feeds) {
        feed->start();
    }
    for(auto feed : unique_output_feeds) {
        feed->start();
    }
    for (auto pivot : pivots) {
        pivot->start();
    }
};
void Pipeline::stop() {
    /*  
        output channel buffers still have residual records
        notify all output feeds that no more input is coming 
        and they should explicitly flush
    */
    for(auto feed : unique_output_feeds) {
        feed->stop();
    }
    for(auto feed : unique_input_feeds) {
        feed->join();
    }
    for(auto feed : unique_output_feeds) {
        feed->join();
    }
};

/*  NucleicAcidAccumulator */

NucleicAcidAccumulator::NucleicAcidAccumulator() :
    count(0),
    min(0),
    max(0),
    sum(0),
    mean(0),
    Q1(0),
    Q3(0),
    IQR(0),
    LW(0),
    RW(0),
    median(0) {
    for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        distribution[i] = 0;
    }
};
inline size_t NucleicAcidAccumulator::quantile(const double portion) {
    size_t position = portion * count;
    size_t cell = 0;
    while (position > 0) {
        if (distribution[cell] >= position) {
            break;
        }
        position -= distribution[cell];
        cell++;
        while (distribution[cell] == 0) {
            cell++;
        }
    }
    return cell;
};
void NucleicAcidAccumulator::finalize() {
    for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        count += distribution[i];
    }
    if (count > 0) {
        for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
            const uint64_t value = distribution[i];
            sum += (value * i);
            if(value != 0) {
                max = i;
                if(min == 0) {
                    min = i;
                }
            }
        }
        mean = double(sum) / double(count);
        median = quantile(0.5);
        Q1 = quantile(0.25);
        Q3 = quantile(0.75);
        IQR = Q3 - Q1;

        double W = Q1 - IQR * 1.5;
        LW = (W < min) ? min : W;

        W = Q3 + IQR * 1.5;
        RW = (W > max) ? max : W;
    }
};
NucleicAcidAccumulator& NucleicAcidAccumulator::operator+=(const NucleicAcidAccumulator& rhs) {
    for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        distribution[i] += rhs.distribution[i];
    }
    return *this;
};

/*  CycleAccumulator */

CycleAccumulator::CycleAccumulator() :
    iupac_nucleic_acid(IUPAC_CODE_SIZE) {
};
void CycleAccumulator::finalize() {
    // accumulate all nucleotide variations in the NO_NUCLEOTIDE accumulative distribution
    for(uint8_t i = 1; i < IUPAC_CODE_SIZE; i++) {
        for(uint8_t p = 0; p < EFFECTIVE_PHRED_RANGE; p++) {
            iupac_nucleic_acid[NO_NUCLEOTIDE].distribution[p] += iupac_nucleic_acid[i].distribution[p];
        }
    }
    for(auto& distribution : iupac_nucleic_acid) {
        distribution.finalize();
    }
};
CycleAccumulator& CycleAccumulator::operator+=(const CycleAccumulator& rhs) {
    for(size_t i = 0; i < iupac_nucleic_acid.size(); i++) {
        iupac_nucleic_acid[i] += rhs.iupac_nucleic_acid[i];
    }
    return *this;
};

/*  SegmentAccumulator */

SegmentAccumulator::SegmentAccumulator() :
    count(0),
    min(0),
    max(0),
    sum(0),
    mean(0) {
    for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        distribution[i] = 0;
    }
};
void SegmentAccumulator::finalize() {
    if(count > 0) {
        mean = sum / double(count);
    }
};
SegmentAccumulator& SegmentAccumulator::operator+=(const SegmentAccumulator& rhs) {
    count += rhs.count;
    sum += rhs.sum;
    min = MIN(min, rhs.min);
    max = MAX(max, rhs.max);
    for(size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        distribution[i] += rhs.distribution[i];
    }
    return *this;
};

/*  FeedAccumulator */

FeedAccumulator::FeedAccumulator() :
    length(0),
    shortest(numeric_limits<uint64_t>::max()) {
    for(size_t i = 0; i < IUPAC_CODE_SIZE; i++) {
        iupac_nucleic_acid_count[i] = 0;
    }
};
FeedAccumulator::~FeedAccumulator() {
    for(auto cycle : cycles) {
        delete cycle;
    }
};
void FeedAccumulator::encode(Document& document, Value& value) const {
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;

    v.SetUint64(shortest);
    value.AddMember("min sequence length", v, allocator);

    v.SetUint64(length);
    value.AddMember("max sequence length", v, allocator);

    Value cycle_quality_report;
    cycle_quality_report.SetObject();

    Value cycle_nucleotide_quality_reports;
    cycle_nucleotide_quality_reports.SetArray();

    for (uint8_t n = 0; n < IUPAC_CODE_SIZE; n++) {
        if(iupac_nucleic_acid_count[n] > 0) {
            Value cycle_quality_distribution;
            cycle_quality_distribution.SetObject();

            Value cycle_count;
            cycle_count.SetArray();

            Value cycle_quality_first_quartile;
            cycle_quality_first_quartile.SetArray();

            Value cycle_quality_third_quartile;
            cycle_quality_third_quartile.SetArray();

            Value cycle_quality_interquartile_range;
            cycle_quality_interquartile_range.SetArray();

            Value cycle_quality_left_whisker;
            cycle_quality_left_whisker.SetArray();

            Value cycle_quality_right_whisker;
            cycle_quality_right_whisker.SetArray();

            Value cycle_quality_min;
            cycle_quality_min.SetArray();

            Value cycle_quality_max;
            cycle_quality_max.SetArray();

            Value cycle_quality_mean;
            cycle_quality_mean.SetArray();

            Value cycle_quality_median;
            cycle_quality_median.SetArray();

            Value cycle_quality_sum;
            cycle_quality_sum.SetArray();

            for (size_t c = 0; c < cycles.size(); c++) {
                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].count);
                cycle_count.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].Q1);
                cycle_quality_first_quartile.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].Q3);
                cycle_quality_third_quartile.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].IQR);
                cycle_quality_interquartile_range.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].LW);
                cycle_quality_left_whisker.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].RW);
                cycle_quality_right_whisker.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].min);
                cycle_quality_min.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].max);
                cycle_quality_max.PushBack(v, allocator);

                v.SetDouble(cycles[c]->iupac_nucleic_acid[n].mean);
                cycle_quality_mean.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].median);
                cycle_quality_median.PushBack(v, allocator);

                v.SetUint64(cycles[c]->iupac_nucleic_acid[n].sum);
                cycle_quality_sum.PushBack(v, allocator);
            }

            cycle_quality_distribution.AddMember("cycle count", cycle_count, allocator);
            cycle_quality_distribution.AddMember("cycle quality first quartile", cycle_quality_first_quartile, allocator);
            cycle_quality_distribution.AddMember("cycle quality third quartile", cycle_quality_third_quartile, allocator);
            cycle_quality_distribution.AddMember("cycle quality interquartile range", cycle_quality_interquartile_range, allocator);
            cycle_quality_distribution.AddMember("cycle quality left whisker", cycle_quality_left_whisker, allocator);
            cycle_quality_distribution.AddMember("cycle quality right whisker", cycle_quality_right_whisker, allocator);
            cycle_quality_distribution.AddMember("cycle quality min", cycle_quality_min, allocator);
            cycle_quality_distribution.AddMember("cycle quality max", cycle_quality_max, allocator);
            cycle_quality_distribution.AddMember("cycle quality mean", cycle_quality_mean, allocator);
            cycle_quality_distribution.AddMember("cycle quality median", cycle_quality_median, allocator);
            cycle_quality_distribution.AddMember("cycle quality sum", cycle_quality_sum, allocator);

            if (n > 0) {
                Value cycle_nucleotide_quality_report;
                cycle_nucleotide_quality_report.SetObject();

                Value nucleotide_count;
                nucleotide_count.SetUint64(iupac_nucleic_acid_count[n]);
                cycle_nucleotide_quality_report.AddMember("nucleotide count", nucleotide_count, allocator);

                Value nucleotide;
                nucleotide.SetString(string(1, BamToAmbiguousAscii[n]).c_str(), allocator);
                cycle_nucleotide_quality_report.AddMember("nucleotide", nucleotide, allocator);

                cycle_nucleotide_quality_report.AddMember("cycle quality distribution", cycle_quality_distribution, allocator);
                cycle_nucleotide_quality_reports.PushBack(cycle_nucleotide_quality_report, allocator);

            } else {
                cycle_quality_report.AddMember("cycle quality distribution", cycle_quality_distribution, allocator);
            }
        }
    }
    value.AddMember("cycle nucleotide quality reports", cycle_nucleotide_quality_reports, allocator);
    value.AddMember("cycle nucleotide quality report", cycle_quality_report, allocator);

    Value average_phred_report;
    average_phred_report.SetObject();

    v.SetDouble(average_phred.min);
    average_phred_report.AddMember("average phred score min", v, allocator);

    v.SetDouble(average_phred.max);
    average_phred_report.AddMember("average phred score max", v, allocator);

    v.SetDouble(average_phred.sum);
    average_phred_report.AddMember("average phred score sum", v, allocator);

    v.SetDouble(average_phred.mean);
    average_phred_report.AddMember("average phred score mean", v, allocator);

    Value SegmentAccumulator;
    SegmentAccumulator.SetArray();

    for (size_t i = 0; i < EFFECTIVE_PHRED_RANGE; i++) {
        v.SetUint64(average_phred.distribution[i]);
        SegmentAccumulator.PushBack(v, allocator);
    }
    average_phred_report.AddMember("average phred score distribution", SegmentAccumulator, allocator);
    value.AddMember("average phred score report", average_phred_report, allocator);
};
void FeedAccumulator::finalize() {
    if(shortest == numeric_limits<uint64_t>::max()) {
        shortest = 0;
    }
    for(auto cycle : cycles) {
        cycle->finalize();
    }
    average_phred.finalize();
};
FeedAccumulator& FeedAccumulator::operator+=(const FeedAccumulator& rhs) {
    if(rhs.length > length) {
        for(size_t i = length; i < rhs.length; i++) {
            cycles.push_back(new CycleAccumulator());
        }
        length = rhs.length;
    }

    shortest = MIN(shortest, rhs.shortest);

    for(size_t code = 0; code < IUPAC_CODE_SIZE; code++) {
        iupac_nucleic_acid_count[code] += rhs.iupac_nucleic_acid_count[code];
    }

    for(size_t i = 0; i < rhs.length; i++) {
        *cycles[i] += *rhs.cycles[i];
    }
    average_phred += rhs.average_phred;
    return *this;
};

/*  PivotAccumulator */

PivotAccumulator::PivotAccumulator(const size_t total_input_segments) :
    count(0),
    filtered_count(0),
    determined_count(0),
    determined_filtered_count(0),
    feed_accumulators(total_input_segments) {
};
void PivotAccumulator::increment(const Pivot& pivot) {
    count++;
    if (pivot.filtered) {
        filtered_count++;
    }
    if (pivot.determined) {
        determined_count++;
    }
    if (pivot.determined && pivot.filtered) {
        determined_filtered_count++;
    }
    for(size_t i = 0; i < feed_accumulators.size(); i++) {
        feed_accumulators[i].increment(pivot.input[i].sequence);
    }
};
void PivotAccumulator::finalize() {
    for(auto& accumulator : feed_accumulators) {
        accumulator.finalize();
    }
};
void PivotAccumulator::encode(Document& document, Value& value, const bool disable_quality_control) const {
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;
    v.SetUint64(count);
    value.AddMember("read count", v, allocator);

    v.SetUint64(filtered_count);
    value.AddMember("filtered read count", v, allocator);

    Value feed_reports;
    feed_reports.SetArray();
    for (auto& accumulator : feed_accumulators) {
        Value feed_report;
        feed_report.SetObject();

        v.SetUint64(count);
        feed_report.AddMember("read count", v, allocator);

        v.SetUint64(filtered_count);
        feed_report.AddMember("filtered read count", v, allocator);

        v.SetDouble(count > 0 ? double(filtered_count) / double(count) : 0.0);
        feed_report.AddMember("filtered read fraction", v, allocator);

        // v.SetString((char*)lane.feeds[i].url.c_str(), allocator);
        // feed_report.AddMember("path", v, allocator);

        if (!disable_quality_control) {
            accumulator.encode(document, feed_report);
        }
        feed_reports.PushBack(feed_report, allocator);
    }
    value.AddMember("fastq quality reports", feed_reports, allocator);
};
PivotAccumulator& PivotAccumulator::operator+=(const PivotAccumulator& rhs) {
    count += rhs.count;
    filtered_count += rhs.filtered_count;
    determined_count += rhs.determined_count;
    determined_filtered_count += rhs.determined_filtered_count;
    for(size_t i = 0; i < feed_accumulators.size(); i++) {
        feed_accumulators[i] += rhs.feed_accumulators[i];
    }
    return *this;
};

/*  ChannelAccumulator */

ChannelAccumulator::ChannelAccumulator(const ChannelSpecification& specification):
    count(0),
    filtered_count(0),
    perfect_count(0),
    yield(0),
    perfect_yield(0),
    filtered_yield(0),
    distance_sum(0),
    probability_sum(0),
    include_filtered(specification.include_filtered),
    feed_accumulators(specification.TC) {
};
void ChannelAccumulator::increment(const Pivot& pivot) {
    count++;
    if (pivot.filtered) {
        filtered_count++;
    }
    if (!pivot.multiplex_distance) {
        perfect_count++;
    }
    if (include_filtered || !pivot.filtered) {
        yield++;
        if (!pivot.multiplex_distance) {
            perfect_yield++;
        }
        if (pivot.filtered) {
            filtered_yield++;
        }
        distance_sum += pivot.multiplex_distance;
        probability_sum += pivot.multiplex_probability;
    }
    for(size_t i = 0; i < feed_accumulators.size(); i++) {
        feed_accumulators[i].increment(pivot.output[i].sequence);
    }
};
void ChannelAccumulator::finalize() {
    for(auto& accumulator : feed_accumulators) {
        accumulator.finalize();
    }
};
void ChannelAccumulator::encode(Document& document, Value& value, const bool disable_quality_control) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    Value v;

    v.SetUint64(count);
    value.AddMember("read count", v, allocator);

    v.SetUint64(filtered_count);
    value.AddMember("filtered read count", v, allocator);

    Value feed_reports;
    feed_reports.SetArray();
    for (auto& accumulator : feed_accumulators) {
        Value feed_report;
        feed_report.SetObject();

        v.SetUint64(count);
        feed_report.AddMember("read count", v, allocator);

        v.SetUint64(perfect_count);
        feed_report.AddMember("perfect read count", v, allocator);

        v.SetUint64(filtered_count);
        feed_report.AddMember("filtered read count", v, allocator);

        v.SetUint64(yield);
        feed_report.AddMember("yield", v, allocator);

        v.SetUint64(perfect_yield);
        feed_report.AddMember("perfect yield", v, allocator);

        v.SetUint64(filtered_yield);
        feed_report.AddMember("filtered yield", v, allocator);

        v.SetDouble(count > 0 ? double(filtered_count) / double(count) : 0.0);
        feed_report.AddMember("filtered read fraction", v, allocator);

        v.SetDouble(count > 0 ? double(yield) / double(count) : 0.0);
        feed_report.AddMember("yield fraction", v, allocator);

        // v.SetDouble(totalYield > 0 ? double(library.yield) / double(totalYield) : 0.0);
        // feed_report.AddMember("effective yield fraction", v, allocator);

        v.SetDouble(yield > 0 ? double(perfect_yield) / double(yield) : 0.0);
        feed_report.AddMember("perfect yield fraction", v, allocator);

        v.SetDouble(yield > 0 ? double(filtered_yield) / double(yield) : 0.0);
        feed_report.AddMember("filtered yield fraction", v, allocator);

        v.SetDouble(probability_sum > 0 ? probability_sum / double(yield) : 0.0);
        feed_report.AddMember("expected confidence", v, allocator);

        v.SetDouble(distance_sum > 0 ? double(distance_sum) / double(yield) : 0.0);
        feed_report.AddMember("expected distance", v, allocator);

        if (!disable_quality_control) {
            accumulator.encode(document, feed_report);
        }
        feed_reports.PushBack(feed_report, allocator);
    }
    value.AddMember("fastq quality reports", feed_reports, allocator);
};
ChannelAccumulator& ChannelAccumulator::operator+=(const ChannelAccumulator& rhs) {
    count += rhs.count;
    filtered_count += rhs.filtered_count;
    perfect_count += rhs.perfect_count;
    yield += rhs.yield;
    perfect_yield += rhs.perfect_yield;
    filtered_yield += rhs.filtered_yield;
    distance_sum += rhs.distance_sum;
    probability_sum += rhs.probability_sum;
    distance_sum += rhs.distance_sum;
    distance_sum += rhs.distance_sum;

    for(size_t i = 0; i < feed_accumulators.size(); i++) {
        feed_accumulators[i] += rhs.feed_accumulators[i];
    }
    return *this;
};

