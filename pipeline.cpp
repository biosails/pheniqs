/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
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

void to_string(const FormatKind& value, string& result) {
    switch(value) {
        case FormatKind::FASTQ: result.assign("FASTQ");      break;
        case FormatKind::HTS:   result.assign("HTS");        break;
        default:                result.assign("UNKNOWN");    break;
    }
};
bool from_string(const char* value, FormatKind& result) {
         if(value == NULL)              result = FormatKind::UNKNOWN;
    else if(!strcmp(value, "FASTQ"))    result = FormatKind::FASTQ;
    else if(!strcmp(value, "HTS"))      result = FormatKind::HTS;
    else                                result = FormatKind::UNKNOWN;

    return (result == FormatKind::UNKNOWN ? false : true);
};
void to_kstring(const FormatKind& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
bool from_string(const string& value, FormatKind& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const FormatKind& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const FormatKind& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< FormatKind >(const Value::Ch* key, FormatKind& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

/*  Pivot */

Pivot::Pivot(Pipeline& pipeline) :
    index(static_cast< int32_t >(pipeline.pivot_array.size())),
    action(pipeline.program_action),
    decoder(pipeline.decoder),
    multiplex_barcode(pipeline.multiplex_segment_cardinality),
    molecular_barcode(pipeline.molecular_segment_cardinality),
    accumulator(pipeline.input_specification),
    pipeline(pipeline),
    disable_quality_control(pipeline.disable_quality_control),
    long_read(pipeline.long_read),
    leading_segment(NULL) {

    /* create input segments */
    input.reserve(pipeline.input_segment_cardinality);
    for(int32_t i = 0; i < pipeline.input_segment_cardinality; i++) {
        input.emplace_back(i, i + 1, pipeline.input_segment_cardinality, pipeline.platform);
    }

    /* create output segments */
    output.reserve(pipeline.output_segment_cardinality);
    for(int32_t i = 0; i < pipeline.output_segment_cardinality; i++) {
        output.emplace_back(i, i + 1, pipeline.output_segment_cardinality, pipeline.platform);
    }

    /* set first output segment READ1 flag ON */
    if(pipeline.output_segment_cardinality > 0) {
        output[0].flag |= uint16_t(HtsFlag::READ1);
    }

    /* set last output segment READ2 flag ON */
    if(pipeline.output_segment_cardinality > 1) {
        output[pipeline.output_segment_cardinality - 1].flag |= uint16_t(HtsFlag::READ2);
    }

    leading_segment = &(input[pipeline.leading_segment_index]);

    /*
        in long read accumulation mode channel statistics
        are collected on the pivots and only at the end summed up
        on the channel. otherwise they are collected directly on the channel
    */
    if(!disable_quality_control && long_read) {
        switch (action) {
            case ProgramAction::DEMULTIPLEX: {
                channel_accumulator_by_barcode.reserve(pipeline.channel_specification_array.size());
                for(const auto& specification : pipeline.channel_specification_array) {
                    if(!specification.multiplex_barcode.empty()) {
                        channel_accumulator_by_barcode.emplace(make_pair(specification.multiplex_barcode, specification));
                    }
                }
                break;
            };
            case ProgramAction::QUALITY: {
                channel_accumulator_by_read_group_id.reserve(pipeline.channel_specification_array.size());
                for(const auto& specification : pipeline.channel_specification_array) {
                    if(specification.rg.ID.l > 0) {
                        channel_accumulator_by_read_group_id.emplace(make_pair(specification.rg, specification));
                    }
                }
                break;
            };
            default:
                break;
        }
    }
    clear();
};
void Pivot::start() {
    pivot_thread = thread(&Pivot::run, this);
};
inline void Pivot::clear() {
    /* Those get assigned either way so no need to reset them
    determined = false;
    filtered = false;
    multiplex_probability = 0;
    conditioned_multiplex_probability = 0;
    multiplex_distance = 0;
    decoded_multiplex_channel = NULL; */
    multiplex_barcode.clear();
    molecular_barcode.clear();
    for(auto& segment : input) segment.clear();
    for(auto& segment : output) segment.clear();
};
void Pivot::run() {
    switch (action) {
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
        default:
            break;
    }
};
inline void Pivot::validate() {
    if(input.size() > 1) {
        /* validate that all segments in the pivot have the same identifier */
        const kstring_t& baseline = input[0].name;
        for(size_t i = 1; i < input.size(); i++) {
            Segment& segment = input[i];
            if((baseline.l != segment.name.l) || strncmp(baseline.s, segment.name.s, baseline.l)) {
                throw SequenceError("read segments out of sync " + string(segment.name.s, segment.name.l) + " and " + string(baseline.s, baseline.l));
            }
        }
    }
};
inline void Pivot::transform() {
    // populate the output
    for(auto& transform : pipeline.template_transform_array ) {
        const Segment& from = input[transform.token.input_segment_index];
        Segment& to = output[transform.output_segment_index];
        to.sequence.append(from.sequence, transform);
    }

    // populate the multiplex barcode
    for(auto& transform : pipeline.multiplex_barcode_transform_array) {
        const Segment& from = input[transform.token.input_segment_index];
        multiplex_barcode.append(transform.output_segment_index, from.sequence, transform);
    }

    // populate the molecular barcode
    for(auto& transform : pipeline.molecular_barcode_transform_array) {
        const Segment& from = input[transform.token.input_segment_index];
        molecular_barcode.append(transform.output_segment_index, from.sequence, transform);
    }

    // assign the pivot qc_fail flag from the leader
    filtered = leading_segment->get_qcfail();

    for(auto& segment : output) {
        ks_put_string(leading_segment->name, segment.name);
        segment.set_qcfail(leading_segment->get_qcfail());
        segment.auxiliary.XI = leading_segment->auxiliary.XI;
    }
};
inline void Pivot::decode_with_pamld() {
    multiplex_distance = pipeline.concatenated_multiplex_barcode_length;
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
    double adjusted(0);
    double compensation(0);
    double sigma(0);
    double y(0);
    double t(0);
    double c(0);
    double p(0);
    int32_t d(0);
    for(auto channel : pipeline.channel_array) {
        channel->multiplex_barcode.accurate_decoding_probability(multiplex_barcode, c, d);
        p = c * channel->concentration;
        y = p - compensation;
        t = sigma + y;
        compensation = (t - sigma) - y;
        sigma = t;
        if(p > adjusted) {
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
    multiplex_probability = adjusted / (sigma + pipeline.adjusted_multiplex_noise_probability);

    /* Check for decoding failure and assign to the undetermined channel if decoding failed */
    if(conditioned_multiplex_probability > pipeline.random_multiplex_barcode_probability &&
        multiplex_probability > pipeline.multiplex_confidence) {
        determined = true;

    } else {
        decoded_multiplex_channel = pipeline.undetermined;
        conditioned_multiplex_probability = 0;
        multiplex_distance = 0;
        multiplex_probability = 1;
    }
};
inline void Pivot::decode_with_mdd() {
    multiplex_distance = pipeline.concatenated_multiplex_barcode_length;
    decoded_multiplex_channel = pipeline.undetermined;
    determined = false;

    /* First try a perfect match to the full barcode sequence */
    auto record = pipeline.channel_by_barcode.find(multiplex_barcode);
    if(record != pipeline.channel_by_barcode.end()) {
        decoded_multiplex_channel = record->second;
        multiplex_distance = 0;
        determined = true;

    } else {
        /* If no exact match was found try error correction */
        for(auto channel : pipeline.channel_array) {
            if(channel->multiplex_barcode.corrected_match(multiplex_barcode, multiplex_distance)) {
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
inline void Pivot::decode_tagged_channel() {
    if(leading_segment->auxiliary.RG.l > 0) {
        string rg(leading_segment->auxiliary.RG.s, leading_segment->auxiliary.RG.l);
        auto record = pipeline.channel_by_read_group_id.find(rg);
        if(record != pipeline.channel_by_read_group_id.end()) {
            decoded_multiplex_channel = record->second;
        }
    }
};
inline void Pivot::encode_auxiliary () {
    if(decoded_multiplex_channel != NULL) {
        for(auto& segment: output) {

            if(decoder == Decoder::PAMLD) {
                segment.auxiliary.DQ = 1.0 - multiplex_probability;
            }

            // sequence auxiliary tags
            segment.sequence.expected_error(segment.auxiliary.EE);

            // pivot auxiliary tags
            segment.auxiliary.set_multiplex_barcode(multiplex_barcode);
            segment.auxiliary.set_molecular_barcode(molecular_barcode);

            // channel auxiliary tags
            ks_put_string(decoded_multiplex_channel->rg.ID, segment.auxiliary.RG);

            /* those don't change between consecutive reads in the same channel
            ks_put_string(decoded_multiplex_channel->rg.LB, &segment.auxiliary.LB);
            ks_put_string(decoded_multiplex_channel->rg.PU, &segment.auxiliary.PU);
            ks_put_string(decoded_multiplex_channel->rg.FS, &segment.auxiliary.FS);
            ks_put_string(decoded_multiplex_channel->rg.PG, &segment.auxiliary.PG);
            ks_put_string(decoded_multiplex_channel->rg.CO, &segment.auxiliary.CO); */
        }
    }
};
inline void Pivot::copy_auxiliary () {
    for(auto& segment: output) {
        segment.auxiliary.imitate(leading_segment->auxiliary);
    }
};
inline void Pivot::push() {
    if(decoded_multiplex_channel != NULL) {
        decoded_multiplex_channel->push(*this);
    }
};
inline void Pivot::increment() {
    if(!disable_quality_control) {
        /* input quality tracking */
        accumulator.increment(filtered, input);

        /* long read track channel quality on the pivot */
        if(long_read) {
            switch (action) {
                case ProgramAction::DEMULTIPLEX: {
                    if(!decoded_multiplex_channel->barcode_key.empty()) {
                        const auto& record = channel_accumulator_by_barcode.find(decoded_multiplex_channel->barcode_key);
                        if(record != channel_accumulator_by_barcode.end()) {
                            record->second.increment(
                                filtered,
                                multiplex_probability,
                                multiplex_distance,
                                output
                            );
                        }
                    }
                    break;
                };
                case ProgramAction::QUALITY: {
                    if(!decoded_multiplex_channel->rg_key.empty()) {
                        const auto& record = channel_accumulator_by_read_group_id.find(decoded_multiplex_channel->rg_key);
                        if(record != channel_accumulator_by_read_group_id.end()) {
                            record->second.increment(
                                filtered,
                                multiplex_probability,
                                multiplex_distance,
                                output
                            );
                        }
                    }
                    break;
                };
                default:
                    break;
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
    barcode_key(specification.multiplex_barcode),
    rg_key(specification.rg),
    concentration(specification.concentration),
    multiplex_barcode(specification.multiplex_barcode),
    disable_quality_control(specification.disable_quality_control),
    long_read(specification.long_read),
    include_filtered(specification.include_filtered),
    undetermined(specification.undetermined),
    writable(!specification.empty()),
    rg(specification.rg),
    channel_accumulator(specification) {
};
Channel::~Channel() {
};
inline void Channel::increment(Pivot& pivot) {
    lock_guard<mutex> channel_lock(channel_mutex);
    channel_accumulator.increment (
        pivot.filtered,
        pivot.multiplex_probability,
        pivot.multiplex_distance,
        pivot.output
    );
};
void Channel::push(Pivot& pivot) {
    if(writable) {
        // acquire a push lock for all feeds in a fixed order
        vector< unique_lock< mutex > > feed_locks;
        feed_locks.reserve(output_feed_by_order.size());
        for(const auto feed : output_feed_by_order) {
            feed_locks.push_back(feed->acquire_push_lock());
        }

        // push the segments to the output feeds
        for(size_t i = 0; i < output_feed_by_segment.size(); i++) {
            output_feed_by_segment[i]->push(pivot.output[i]);
        }

        // release the locks on the feeds in reverse order
        for(auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
            feed_lock->unlock();
        }
    }

    // for short reads the channels track their own quality
    if(!disable_quality_control && !long_read) {
        increment(pivot);
    }
};
void Channel::encode(Document& document, Value& value) const {
    channel_accumulator.encode(document, value);
};
void Channel::finalize(const PipelineAccumulator& pipeline_accumulator) {
    channel_accumulator.finalize(pipeline_accumulator);
};

/*  Pipeline */

Pipeline::Pipeline(Environment& environment) :
    environment(environment),
    instruction(environment.instruction()),
    program_action(environment.program_action()),
    decoder(decode_value_by_key< Decoder >("decoder", instruction)),
    platform(decode_value_by_key< Platform >("platform", instruction)),
    leading_segment_index(decode_value_by_key< int32_t >("leading segment index", instruction)),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", instruction)),
    long_read(decode_value_by_key< bool >("long read", instruction)),
    include_filtered(decode_value_by_key< bool >("include filtered", instruction)),
    adjusted_multiplex_noise_probability(decode_value_by_key< double >("adjusted multiplex noise probability", instruction)),
    random_multiplex_barcode_probability(decode_value_by_key< double >("random multiplex barcode probability", instruction)),
    multiplex_confidence(decode_value_by_key< double >("multiplex confidence", instruction)),
    threads(decode_value_by_key< int32_t >("threads", instruction)),
    input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", instruction)),
    output_segment_cardinality(decode_value_by_key< int32_t >("output segment cardinality", instruction)),
    multiplex_segment_cardinality(decode_value_by_key< int32_t >("multiplex segment cardinality", instruction)),
    molecular_segment_cardinality(decode_value_by_key< int32_t >("molecular segment cardinality", instruction)),
    concatenated_multiplex_barcode_length(decode_value_by_key< int32_t >("concatenated multiplex barcode length", instruction)),
    concatenated_molecular_barcode_length(decode_value_by_key< int32_t >("concatenated molecular barcode length", instruction)),
    end_of_input(false),
    thread_pool({NULL, 0}),
    undetermined(NULL),
    input_accumulator(NULL),
    output_accumulator(NULL) {

    decode_value_by_key< vector< Token > >("token", token_array, instruction);
    decode_transform_array_by_key("template", template_transform_array, instruction, token_array);
    decode_transform_array_by_key("multiplex barcode", multiplex_barcode_transform_array, instruction, token_array);
    decode_transform_array_by_key("molecular barcode", molecular_barcode_transform_array, instruction, token_array);

    thread_pool.pool = hts_tpool_init(threads);
    if(!thread_pool.pool) {
        throw InternalError("error creating thread pool");
    }
};
Pipeline::~Pipeline() {
    hts_tpool_destroy(thread_pool.pool);
    for(auto feed : input_feed_by_index) {
        delete feed;
    }
    for(auto feed : output_feed_by_index) {
        delete feed;
    }
    input_feed_by_segment.clear();
    input_feed_by_index.clear();
    output_feed_by_index.clear();

    for(auto pivot : pivot_array) {
        delete pivot;
    }
    for(auto channel : channel_array) {
        delete channel;
    }
    delete undetermined;
    undetermined = NULL;
    channel_array.clear();
    channel_by_barcode.clear();
    channel_by_read_group_id.clear();

    delete input_accumulator;
    delete output_accumulator;
};
void Pipeline::execute() {
    switch (program_action) {
        case ProgramAction::DEMULTIPLEX: {
            validate_url_accessibility();
            load_input();
            load_output();
            break;
        };
        case ProgramAction::QUALITY: {
            // probe(environment.input_url);
            // environment.load_transformation();
            // environment.load_channels();
            // environment.load_input_specification();
            load_output();
            break;
        };
        default:
            break;
    }
    load_pivot_array();
    start();
    for(auto pivot : pivot_array) {
        pivot->pivot_thread.join();
    }
    stop();
    if(!disable_quality_control) {
        finalize();
        encode_report(cout);
    }
};
bool Pipeline::pull(Pivot& pivot) {
    vector< unique_lock< mutex > > feed_locks;
    feed_locks.reserve(input_feed_by_index.size());

    // acquire a pull lock for all feeds in a fixed order
    for(const auto feed : input_feed_by_index) {
        feed_locks.push_back(feed->acquire_pull_lock());
    }

    // pull into pivot input segments from input feeds
    for(size_t i = 0; i < pivot.input.size(); i++) {
        if(!input_feed_by_segment[i]->pull(pivot.input[i])) {
            end_of_input = true;
        }
    }

    // release the locks on the feeds in reverse order
    for(auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }
    return !end_of_input;
};
void Pipeline::validate_url_accessibility() {
    URL url;
    Value::ConstMemberIterator collection;

    collection = instruction.FindMember("input feed");
    if(collection != instruction.MemberEnd()) {
        for(auto& element : collection->value.GetArray()) {
            if(decode_file_url_by_key("url", url, IoDirection::IN, element)) {
                if(!url.is_readable()) {
                    throw IOError("could not open " + string(url) + " for reading");
                }
            }
        }
    }

    collection = instruction.FindMember("output feed");
    if(collection != instruction.MemberEnd()) {
        for(auto& element : collection->value.GetArray()) {
            if(decode_file_url_by_key("url", url, IoDirection::IN, element)) {
                if(!url.is_writable()) {
                    throw IOError("could not open " + string(url) + " for writing");
                }
            }
        }
    }
};
void Pipeline::load_input() {
    /* decode the enumerated feed specification array */
    list< FeedSpecification* > feed_specification_array;
    decode_value_by_key< list< FeedSpecification* > >("input feed", feed_specification_array, instruction);

    for(auto specification : feed_specification_array) {
        specification->probe();
    };

    /* load a map of feed specification by URL */
    unordered_map< URL, FeedSpecification* > feed_specification_by_url;
    for(auto specification : feed_specification_array) {
        feed_specification_by_url.emplace(make_pair(specification->url, specification));
        // feed_specification_by_url[specification->url] = specification;
    };

    /*  load the input specification */
    decode_value_by_key< InputSpecification >(NULL, input_specification, instruction);
    input_specification.feed_specification_by_segment.reserve(input_specification.url_by_segment.size());

    /*  populate the feed_specification_by_segment in the input_specification */
    for(const auto& url : input_specification.url_by_segment) {
        const auto& record = feed_specification_by_url.find(url);
        if(record != feed_specification_by_url.end()) {
            input_specification.feed_specification_by_segment.push_back(record->second);
        } else {
            throw InternalError("missing feed specification for URL " + string(url) + " referenced in input specification segment array");
        }
    }

    /*  load input feeds from specification and populate the input_feed_by_url and input_feed_by_index arrays */
    unordered_map< URL, Feed* > input_feed_by_url(feed_specification_array.size());
    for(auto specification : feed_specification_array) {
        Feed* feed(NULL);
        FormatKind kind(format_kind_from_type(specification->url.type()));
        switch(kind) {
            case FormatKind::FASTQ: {
                feed = new FastqFeed(*specification);
                break;
            };
            case FormatKind::HTS: {
                feed = new HtsFeed(*specification);
                break;
            };
            default: {
                throw InternalError("unknown input format " + string(specification->url));
                break;
            };
        }
        feed->set_thread_pool(&thread_pool);
        input_feed_by_index.push_back(feed);
        // input_feed_by_url[specification->url] = feed;
        input_feed_by_url.emplace(make_pair(specification->url, feed));
    }

    /*  populate the input_feed_by_segment array */
    input_feed_by_segment.reserve(input_specification.url_by_segment.size());
    for(const auto& url : input_specification.url_by_segment) {
        const auto& record = input_feed_by_url.find(url);
        if(record != input_feed_by_url.end()) {
            input_feed_by_segment.push_back(record->second);
        } else {
            throw InternalError("missing feed for URL " + string(url) + " referenced in input specification segment array");
        }
    }

    if(!disable_quality_control) {
        input_accumulator = new PivotAccumulator(input_specification);
    }
};
void Pipeline::load_output() {
    /* decode the enumerated feed specification array */
    list< FeedSpecification* > feed_specification_array;
    decode_value_by_key< list< FeedSpecification* > >("output feed", feed_specification_array, instruction);

    for(auto specification : feed_specification_array) {
        specification->probe();
    };

    /* load a map of feed specification by URL */
    unordered_map< URL, FeedSpecification* > feed_specification_by_url;
    for(auto specification : feed_specification_array) {
        feed_specification_by_url.emplace(make_pair(specification->url, specification));
        // feed_specification_by_url[specification->url] = specification;
    };

    /*  load the channel specification array */
    decode_value_by_key< list< ChannelSpecification > >("channel", channel_specification_array, instruction);

    /*  populate the feed_specification_by_segment in each channel specification */
    for(auto& specification : channel_specification_array) {
        specification.feed_specification_by_segment.reserve(specification.url_by_segment.size());
        for(const auto& url : specification.url_by_segment) {
            const auto& record = feed_specification_by_url.find(url);
            if(record != feed_specification_by_url.end()) {
                specification.feed_specification_by_segment.push_back(record->second);
            } else {
                throw InternalError("missing feed specification for URL " + string(url) + " referenced in output specification segment array");
            }
        }
    }

    /*  load output feeds from specification and populate the output_feed_by_url and output_feed_by_index arrays */
    unordered_map< URL, Feed* > output_feed_by_url(feed_specification_array.size());
    for(auto specification : feed_specification_array) {
        Feed* feed(NULL);
        FormatKind kind(format_kind_from_type(specification->url.type()));
        switch(kind) {
            case FormatKind::FASTQ: {
                feed = new FastqFeed(*specification);
                break;
            };
            case FormatKind::HTS: {
                feed = new HtsFeed(*specification);
                break;
            };
            default: {
                throw InternalError("unknown input format " + string(specification->url));
                break;
            };
        }
        feed->set_thread_pool(&thread_pool);
        output_feed_by_index.push_back(feed);
        output_feed_by_url.emplace(make_pair(specification->url, feed));
        // output_feed_by_url[specification->url] = feed;
    }

    /* load channel channel_array */
    for(const auto& specification : channel_specification_array) {
        Channel* channel = new Channel(*this, specification);
        /*  the undetermined channel is not present in the channels array
            and is referenced separately by the undetermined pointer */
        if(!specification.undetermined) {
            channel_array.emplace_back(channel);
        } else {
            undetermined = channel;
        }

        set< URL > unique_feed_url;
        channel->output_feed_by_segment.reserve(specification.url_by_segment.size());
        for(const auto& url : specification.url_by_segment) {
            const auto& record = output_feed_by_url.find(url);
            if(record != output_feed_by_url.end()) {
                channel->output_feed_by_segment.push_back(record->second);
                if(!unique_feed_url.count(url)) {
                    unique_feed_url.emplace(url);
                    channel->output_feed_by_order.push_back(record->second);
                }
            } else {
                throw InternalError("missing feed for URL " + string(url) + " referenced in output specification segment array");
            }
        }

        /* channel by read group id lookup table */
        if(!channel->rg_key.empty()) {
            if(channel_by_read_group_id.find(channel->rg_key) == channel_by_read_group_id.end()) {
                channel_by_read_group_id.emplace(make_pair(channel->rg_key, channel));
                // channel_by_read_group_id[channel->rg_key] = channel;
            }
        }

        /* channel by concatenated barcode string lookup table */
        if(!channel->barcode_key.empty()) {
            if(channel_by_barcode.find(channel->barcode_key) == channel_by_barcode.end()) {
                channel_by_read_group_id.emplace(make_pair(channel->barcode_key, channel));
                // channel_by_barcode[channel->barcode_key] = channel;
            }
        }
    }

    if(!disable_quality_control) {
        output_accumulator = new PipelineAccumulator(decoder);
    }
};
void Pipeline::load_pivot_array() {
    pivot_array.reserve(threads);
    for(int i = 0; i < threads; i++) {
        pivot_array.push_back(new Pivot(*this));
    }
};
void Pipeline::start() {
    for(auto feed : input_feed_by_index) {
        feed->open();
    }
    for(auto feed : output_feed_by_index) {
        feed->open();
    }
    for(auto feed : input_feed_by_index) {
        feed->start();
    }
    for(auto feed : output_feed_by_index) {
        feed->start();
    }
    for(auto pivot : pivot_array) {
        pivot->start();
    }
};
void Pipeline::stop() {
    /*
        output channel buffers still have residual records
        notify all output feeds that no more input is coming
        and they should explicitly flush
    */
    for(auto feed : output_feed_by_index) {
        feed->stop();
    }
    for(auto feed : input_feed_by_index) {
        feed->join();
    }
    for(auto feed : output_feed_by_index) {
        feed->join();
    }
};
void Pipeline::finalize() {
    /* collect statistics from all parallel pivots */
    for(const auto pivot : pivot_array) {
        *input_accumulator += pivot->accumulator;

        if(long_read) {
            /* collect output statistics from pivots on the channel accumulators */
            switch (program_action) {
                case ProgramAction::DEMULTIPLEX: {
                    for(const auto& record : pivot->channel_accumulator_by_barcode) {
                        const auto& channel = channel_by_barcode.find(record.first);
                        if(channel != channel_by_barcode.end()) {
                            channel->second->channel_accumulator += record.second;
                        }
                    }
                    break;
                };
                case ProgramAction::QUALITY: {
                    for(const auto& record : pivot->channel_accumulator_by_read_group_id) {
                        const auto& channel = channel_by_read_group_id.find(record.first);
                        if(channel != channel_by_read_group_id.end()) {
                            channel->second->channel_accumulator += record.second;
                        }
                    }
                    break;
                };
                default:
                    break;
            }
        }
    }
    input_accumulator->finalize();

    /* collect statistics from all output channels */
    for(auto channel : channel_array) {
        output_accumulator->collect(channel->channel_accumulator);
    }
    if(undetermined != NULL) {
        output_accumulator->collect(undetermined->channel_accumulator);
    }

    /* finalize the pipeline accumulator */
    output_accumulator->finalize();

    /* finalize the channel accumulators, now that we know the totals */
    if(undetermined != NULL) {
        undetermined->finalize(*output_accumulator);
    }
    for(auto channel : channel_array) {
        channel->finalize(*output_accumulator);
    }
};
void Pipeline::encode_report(ostream& o) const {
    switch (program_action) {
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
    o << endl;
};
void Pipeline::encode_quality_report(ostream& o) const {
    Document document(kObjectType);
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;

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
void Pipeline::encode_demultiplex_report(ostream& o) const {
    Document document(kObjectType);
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;

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
void Pipeline::encode_input_report(Document& document, Value& value) const {
    value.SetObject();
    input_accumulator->encode(document, value);
};
void Pipeline::encode_output_report(Document& document, Value& value) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    value.SetObject();

    output_accumulator->encode(document, value);

    Value channel_reports(kArrayType);
    /* encode the undetermined channel first */
    if(undetermined != NULL) {
        Value channel_report(kObjectType);
        undetermined->encode(document, channel_report);
        channel_reports.PushBack(channel_report, allocator);
    }

    /* encode the classified channels */
    for(auto channel : channel_array) {
        Value channel_report(kObjectType);
        channel->encode(document, channel_report);
        channel_reports.PushBack(channel_report, allocator);
    }
    value.AddMember("read group quality reports", channel_reports.Move(), allocator);
};
void Pipeline::probe(const URL& url) {
    // if(!url.empty()) {
    //     FeedSpecification* specification = environment.discover_feed(url, IoDirection::IN);
    //     specification->set_resolution(1);
    //     specification->set_capacity(environment.buffer_capacity);
    //     specification->probe();

    //     Feed* feed = load_feed(specification);
    //     feed->open();
    //     FormatKind kind(format_kind_from_type(specification->url.type()));
    //     switch(kind) {
    //         case FormatKind::FASTQ: {
    //             feed->replenish();
    //             Segment segment(specification->platform);
    //             uint64_t count = 0;
    //             if(feed->peek(segment, count)) {
    //                 count++;
    //                 string primary(segment.name.s, segment.name.l);
    //                 while(feed->peek(segment, count) && primary.size() == segment.name.l && strncmp(primary.c_str(), segment.name.s, segment.name.l)) {
    //                     count++;
    //                 }
    //             }
    //             environment.total_output_segments = count;
    //             environment.total_input_segments = count;
    //             break;
    //         };
    //         case FormatKind::HTS:{
    //             feed->replenish();
    //             Segment segment(specification->platform);
    //             feed->peek(segment, 0);
    //             environment.total_input_segments = segment.get_total_segments();
    //             environment.total_output_segments = segment.get_total_segments();

    //             const HtsHeader& header = ((HtsFeed*)feed)->get_header();
    //             for(const auto& record : header.read_group_by_id) {
    //                 ChannelSpecification* channel = environment.load_channel_from_rg(record.second);
    //                 channel->TC = environment.total_output_segments;
    //             }
    //             break;
    //         };
    //         default: break;
    //     }
    //     specification->set_resolution(environment.total_input_segments);
    //     specification->set_capacity(environment.buffer_capacity * specification->resolution);
    //     feed->calibrate(specification);
    //     environment.calibrate(url);
    // }
};
void Pipeline::calibrate(const URL& url) {
    // _input_specification.decoder = _decoder;
    // _input_specification.disable_quality_control = _disable_quality_control;
    // _input_specification.long_read = _long_read;
    // _input_specification.include_filtered = _include_filtered;

    // _input_specification.input_urls.reserve(total_input_segments);
    // for(size_t i = 0; i < total_input_segments; i++) {
    //     _input_specification.input_urls.push_back(url);
    //     token_patterns.emplace_back(to_string(i) + "::");
    //     template_patterns.emplace_back(to_string(i));
    // }
};
ChannelSpecification* Pipeline::load_channel_from_rg(const HeadRGAtom& rg) {
    // ChannelSpecification* specification = new ChannelSpecification(channel_specifications.size());
    // specification->rg = rg;
    // channel_specifications.push_back(specification);
    // return specification;
    return NULL;
};
