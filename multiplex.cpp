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

#include "multiplex.h"

static int32_t compute_inheritence_depth(const string& key, const unordered_map< string, Value* >& object_by_key, Document& document) {
    int32_t depth(0);
    auto record = object_by_key.find(key);
    if(record != object_by_key.end()) {
        Value* value = record->second;
        if(!decode_value_by_key("depth", depth, *value)) {
            string base_key;
            if(decode_value_by_key("base", base_key, *value)) {
                if(base_key != key) {
                    depth = compute_inheritence_depth(base_key, object_by_key, document) + 1;
                    encode_key_value("depth", depth, *value, document);
                } else { throw ConfigurationError(key + " references itself as parent"); }
            } else {
                encode_key_value("depth", depth, *value, document);
            }
        }
    } else { throw ConfigurationError("referencing an unknown parent " + key); }
    return depth;
};

MultiplexJob::MultiplexJob(Document& operation) try :
    Job(operation),
    decoder_repository_query("/decoder"),
    end_of_input(false),
    thread_pool({NULL, 0}) {

    } catch(Error& error) {
        error.push("MultiplexJob");
        throw;
};
MultiplexJob::~MultiplexJob() {
    if(thread_pool.pool != NULL) {
        hts_tpool_destroy(thread_pool.pool);
    }
    for(auto feed : input_feed_by_index) {
        delete feed;
    }
    for(auto feed : output_feed_by_index) {
        delete feed;
    }
    input_feed_by_segment.clear();
    input_feed_by_index.clear();
    output_feed_by_index.clear();
};
void MultiplexJob::assemble() {
    Job::assemble();
    apply_inheritance();
    clean();
};
void MultiplexJob::compile() {
    /* overlay on top of the default configuration */
    apply_default();

    /* overlay interactive parameters on top of the configuration */
    apply_interactive();

    /* compile a PG SAM header ontology with details about pheniqs */
    compile_PG();

    compile_input();
    compile_barcode_decoding();
    compile_output();

    /* Remove the decoder repository, it is no longer needed in the compiled instruction */
    ontology.RemoveMember("decoder");

    Job::compile();
};
void MultiplexJob::validate() {
    Job::validate();

    uint8_t input_phred_offset;
    if(decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, ontology)) {
        if(input_phred_offset > MAX_PHRED_VALUE || input_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("input phred offset out of range " + to_string(input_phred_offset));
        }
    }

    uint8_t output_phred_offset;
    if(decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, ontology)) {
        if(output_phred_offset > MAX_PHRED_VALUE || output_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("output phred offset out of range " + to_string(output_phred_offset));
        }
    }

    validate_decoder_group("multiplex");
    validate_decoder_group("molecular");
    validate_decoder_group("cellular");
};
void MultiplexJob::load() {
    validate_url_accessibility();
    load_thread_pool();
    load_input();
    load_output();
    load_pivot();
};
void MultiplexJob::execute() {
    load();
    start();
    stop();
    finalize();
};

void MultiplexJob::start() {
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
    for(auto& pivot : pivot_array) {
        pivot.start();
    }
    for(auto& pivot : pivot_array) {
        pivot.join();
    }
};
void MultiplexJob::stop() {
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
void MultiplexJob::finalize() {
    Value value;
    value.CopyFrom(ontology, report.GetAllocator());
    report.AddMember(Value("job", report.GetAllocator()).Move(), value.Move(), report.GetAllocator());

    InputAccumulator input_accumulator(ontology);
    OutputAccumulator output_accumulator(find_value_by_key("multiplex", ontology));
    for(auto& pivot : pivot_array) {
        input_accumulator += pivot.input_accumulator;
        output_accumulator += pivot.output_accumulator;
    }
    input_accumulator.finalize();
    output_accumulator.finalize();
    encode_key_value("demultiplex output report", output_accumulator, report, report);
    encode_key_value("demultiplex input report", input_accumulator, report, report);

    clean_json_value(report, report);
    sort_json_value(report, report);
};
bool MultiplexJob::pull(Read& read) {
    vector< unique_lock< mutex > > feed_locks;
    feed_locks.reserve(input_feed_by_index.size());

    /* acquire a pull lock for all input feeds in a fixed order */
    for(const auto feed : input_feed_by_index) {
        feed_locks.push_back(feed->acquire_pull_lock());
    }

    /* pull into pivot input segments from input feeds */
    for(size_t i(0); i < read.segment_cardinality(); ++i) {
        if(!input_feed_by_segment[i]->pull(read[i])) {
            end_of_input = true;
        }
    }

    /* release the locks on the input feeds in reverse order */
    for(auto feed_lock(feed_locks.rbegin()); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }
    return !end_of_input;
};
void MultiplexJob::describe(ostream& o) const {
    print_global_instruction(o);
    print_input_instruction(o);
    print_transform_instruction(o);
    print_multiplex_instruction(o);
    print_molecular_instruction(o);
    print_cellular_instruction(o);
};

void MultiplexJob::apply_inheritance() {
    apply_repository_inheritence("decoder", ontology, ontology);
    apply_topic_inheritance("multiplex");
    apply_topic_inheritance("molecular");
    apply_topic_inheritance("cellular");
};
void MultiplexJob::apply_repository_inheritence(const Value::Ch* key, Value& container, Document& document) {
    if(container.IsObject()) {
        Value::MemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(reference->value.IsObject()) {
                /*  map each object by key */
                unordered_map< string, Value* > object_by_key;
                for(auto& record : reference->value.GetObject()) {
                    if(!record.value.IsNull()) {
                        record.value.RemoveMember("depth");
                        object_by_key.emplace(make_pair(string(record.name.GetString(), record.name.GetStringLength()), &record.value));
                    }
                }

                /* compute the inheritence depth of each object and keep track of the max depth */
                int32_t max_depth(0);
                for(auto& record : object_by_key) {
                    try {
                        max_depth = max(max_depth, compute_inheritence_depth(record.first, object_by_key, document));
                    } catch(ConfigurationError& error) {
                        throw CommandLineError(record.first + " is " + error.message);
                    }
                }

                /* apply object inheritence down the tree */
                int32_t depth(0);
                for(int32_t i(1); i <= max_depth; ++i) {
                    for(auto& record : object_by_key) {
                        Value* value = record.second;
                        if(decode_value_by_key("depth", depth, *value) && depth == i) {
                            string base;
                            if(decode_value_by_key< string >("base", base, *value)) {
                                merge_json_value(*object_by_key[base], *value, document);
                                value->RemoveMember("base");
                            }
                        }
                    }
                }

                for(auto& record : object_by_key) {
                    record.second->RemoveMember("depth");
                }
            }
        }
    }
};
void MultiplexJob::apply_topic_inheritance(const Value::Ch* key) {
    if(ontology.IsObject()) {
        Value::MemberIterator reference = ontology.FindMember(key);
        if(reference != ontology.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    try {
                        apply_decoder_inheritance(reference->value);
                    } catch(ConfigurationError& error) {
                        error.message.insert(0, string(key) + " decoder : ");
                        throw error;
                    }
                } else if(reference->value.IsArray()) {
                    int32_t index(0);
                    try {
                        for(auto& element : reference->value.GetArray()) {
                            apply_decoder_inheritance(element);
                            ++index;
                        }
                    } catch(ConfigurationError& error) {
                        error.message.insert(0, string(key) + " decoder at " + to_string(index) + " : ");
                        throw error;
                    }
                }
            }
        }
    }
};
void MultiplexJob::apply_decoder_inheritance(Value& value) {
    if(value.IsObject()) {
        string base;
        Value::ConstMemberIterator reference;
        if(decode_value_by_key< string >("base", base, value)) {
            const Value* decoder_repository(decoder_repository_query.Get(ontology));
            if(decoder_repository != NULL) {
                reference = decoder_repository->FindMember(base.c_str());
                if(reference != decoder_repository->MemberEnd()) {
                    merge_json_value(reference->value, value, ontology);
                } else { throw ConfigurationError("reference to an unknown base " + base); }
            }
        }
        value.RemoveMember("base");
        clean_json_value(value, ontology);
    }
};

void MultiplexJob::compile_PG() {
    Value PG(kObjectType);

    string buffer;
    if(decode_value_by_key< string >("application name", buffer, ontology)) {
        encode_key_value("ID", buffer, PG, ontology);
    }
    if(decode_value_by_key< string >("application name", buffer, ontology)) {
        encode_key_value("PN", buffer, PG, ontology);
    }
    if(decode_value_by_key< string >("full command", buffer, ontology)) {
        encode_key_value("CL", buffer, PG, ontology);
    }
    if(decode_value_by_key< string >("previous application", buffer, ontology)) {
        encode_key_value("PP", buffer, PG, ontology);
    }
    if(decode_value_by_key< string >("application description", buffer, ontology)) {
        encode_key_value("DS", buffer, PG, ontology);
    }
    if(decode_value_by_key< string >("application version", buffer, ontology)) {
        encode_key_value("VN", buffer, PG, ontology);
    }
    ontology.AddMember(Value("program", ontology.GetAllocator()).Move(), PG.Move(), ontology.GetAllocator());
};
void MultiplexJob::compile_input() {
    /* Populate the input_feed_by_index and input_feed_by_segment arrays */

    Platform platform(decode_value_by_key< Platform >("platform", ontology));
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    uint8_t input_phred_offset(decode_value_by_key< uint8_t >("input phred offset", ontology));

    expand_url_value_by_key("base input url", ontology, ontology);
    URL base(decode_value_by_key< URL >("base input url", ontology));

    expand_url_array_by_key("input", ontology, ontology, IoDirection::IN);
    relocate_url_array_by_key("input", ontology, ontology, base);
    list< URL > explicit_url_array(decode_value_by_key< list< URL > >("input", ontology));

    int32_t feed_index(0);
    int32_t input_segment_cardinality(0);
    if(sense_input_layout()) {
        /*  create a proxy for each unique feed
            feed_proxy_by_index is ordered by the first appearance of the url */
        list< FeedProxy > feed_proxy_by_index;
        unordered_map< URL, FeedProxy* > feed_proxy_by_url;
        for(const auto& url : explicit_url_array) {
            if(feed_proxy_by_url.count(url) == 0) {
                Value element(kObjectType);
                encode_key_value("index", feed_index, element, ontology);
                encode_key_value("url", url, element, ontology);
                encode_key_value("direction", IoDirection::IN, element, ontology);
                encode_key_value("platform", platform, element, ontology);
                encode_key_value("capacity", buffer_capacity, element, ontology);
                encode_key_value("resolution", 1, element, ontology);
                encode_key_value("phred offset", input_phred_offset, element, ontology);
                feed_proxy_by_index.emplace_back(element);
                feed_proxy_by_url.emplace(make_pair(url, &feed_proxy_by_index.back()));
                ++feed_index;
            }
        }

        /*  Detect the resolution of every feed
            Count the number or consecutive reads with identical read id
            The resolution is the number of segments of each read interleaved into the file
            the read id of all interleaved blocks from all input feeds must match.
            The sum of all uniqe input feed resolutions is the input segment cardinality */
        unordered_map< URL, Feed* > feed_by_url;
        unordered_map< URL, string > read_id_by_url;
        for(auto& proxy : feed_proxy_by_index) {
            Feed* feed(NULL);
            int32_t resolution(0);
            proxy.probe();
            Segment segment;
            string feed_read_id;
            switch(proxy.kind()) {
                case FormatKind::FASTQ: {
                    FastqFeed* fastq_feed = new FastqFeed(proxy);
                    fastq_feed->set_thread_pool(&thread_pool);
                    fastq_feed->open();
                    fastq_feed->replenish();
                    if(fastq_feed->peek(segment, resolution)) {
                        ++resolution;
                        feed_read_id.assign(segment.name.s, segment.name.l);
                        while (
                            fastq_feed->peek(segment, resolution) &&
                            feed_read_id.size() == segment.name.l &&
                            strncmp(feed_read_id.c_str(), segment.name.s, segment.name.l)
                        ) { ++resolution; }
                    }
                    fastq_feed->calibrate_resolution(resolution);
                    feed = fastq_feed;
                    break;
                };
                case FormatKind::HTS: {
                    HtsFeed* hts_feed = new HtsFeed(proxy);
                    hts_feed->set_thread_pool(&thread_pool);
                    hts_feed->open();
                    hts_feed->replenish();
                    if(hts_feed->peek(segment, resolution)) {
                        feed_read_id.assign(segment.name.s, segment.name.l);
                        resolution = segment.total_segments();
                    }
                    hts_feed->calibrate_resolution(resolution);
                    /*
                        const HtsHeader& header = ((HtsFeed*)feed)->get_header();
                        for(const auto& record : header.read_group_by_id) {}
                    */
                    feed = hts_feed;
                    break;
                };
                case FormatKind::DEV_NULL: {
                    throw ConfigurationError("/dev/null can not be used for input");
                    break;
                };
                default: {
                    throw ConfigurationError("unknown input format " + string(proxy.url));
                    break;
                };
            }

            read_id_by_url[proxy.url] = feed_read_id;
            input_feed_by_index.push_back(feed);
            feed_by_url.emplace(make_pair(proxy.url, feed));
            input_segment_cardinality += resolution;
        };

        /* Check that the read id of all interleaved blocks from all input feeds match. */
        if(input_segment_cardinality > 1) {
            string anchor_read_id;
            URL anchor_url;
            for(auto& record : read_id_by_url) {
                if(anchor_read_id.empty()) {
                    anchor_url = record.first;
                    anchor_read_id = record.second;
                } else if(anchor_read_id != record.second) {
                    throw ConfigurationError(string(anchor_url) + " and " + string(record.second) + " are out of sync");
                }
            }
        }
        encode_key_value("input segment cardinality", input_segment_cardinality, ontology, ontology);

        input_feed_by_segment.reserve(input_segment_cardinality);
        for(auto& feed : input_feed_by_index) {
            for(int32_t i(0); i < feed->resolution(); ++i) {
                input_feed_by_segment.push_back(feed);
            }
        }

        list< URL > feed_url_by_segment;
        for(auto& feed : input_feed_by_segment) {
            feed_url_by_segment.push_back(feed->url);
        }
        encode_key_value("input", feed_url_by_segment, ontology, ontology);
        encode_key_value("input feed", input_feed_by_index, ontology, ontology);
        encode_key_value("input feed by segment", input_feed_by_segment, ontology, ontology);

    } else {
        /* encode the input segment cardinality */
        input_segment_cardinality = static_cast< int32_t>(explicit_url_array.size());
        encode_key_value("input segment cardinality", input_segment_cardinality, ontology, ontology);

        list< URL > feed_url_by_index;
        map< URL, int > feed_resolution;
        for(const auto& url : explicit_url_array) {
            auto record = feed_resolution.find(url);
            if(record == feed_resolution.end()) {
                feed_resolution[url] = 1;
                feed_url_by_index.push_back(url);
            } else { ++(record->second); }
        }

        unordered_map< URL, Value > feed_ontology_by_url;
        for(const auto& url : feed_url_by_index) {
            Value element(kObjectType);
            int resolution(feed_resolution[url]);
            encode_key_value("index", feed_index, element, ontology);
            encode_key_value("url", url, element, ontology);
            encode_key_value("direction", IoDirection::IN, element, ontology);
            encode_key_value("platform", platform, element, ontology);
            encode_key_value("capacity", buffer_capacity, element, ontology);
            encode_key_value("resolution", resolution, element, ontology);
            encode_key_value("phred offset", input_phred_offset, element, ontology);
            feed_ontology_by_url.emplace(make_pair(url, move(element)));
            ++feed_index;
        }

        Value feed_by_segment_array(kArrayType);
        for(const auto& url : explicit_url_array) {
            Value feed_ontology(feed_ontology_by_url[url], ontology.GetAllocator());
            feed_by_segment_array.PushBack(feed_ontology.Move(), ontology.GetAllocator());
        }
        ontology.RemoveMember("input feed by segment");
        ontology.AddMember("input feed by segment", feed_by_segment_array.Move(), ontology.GetAllocator());

        Value feed_by_index_array(kArrayType);
        for(const auto& url : feed_url_by_index) {
            feed_by_index_array.PushBack(feed_ontology_by_url[url].Move(), ontology.GetAllocator());
        }
        feed_ontology_by_url.clear();

        ontology.RemoveMember("input feed");
        ontology.AddMember("input feed", feed_by_index_array.Move(), ontology.GetAllocator());
    }

    /*  validate leading_segment_index */
    int32_t leading_segment_index(decode_value_by_key< int32_t >("leading segment index", leading_segment_index, ontology));
    if(leading_segment_index >= input_segment_cardinality) {
        throw ConfigurationError("leading segment index " + to_string(leading_segment_index) + " references non existing input segment");
    }
};
void MultiplexJob::compile_barcode_decoding() {
    compile_topic("multiplex");
    compile_topic("molecular");
    compile_topic("cellular");
};
void MultiplexJob::compile_topic(const Value::Ch* key) {
    Value::MemberIterator reference = ontology.FindMember(key);
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            /* aseemble a default decoder */
            string decoder_projection_uri(key);
            decoder_projection_uri.append(":decoder");
            const Value* decoder_projection(find_projection(decoder_projection_uri));

            Value decoder_template(kObjectType);
            if(decoder_projection != NULL && !decoder_projection->IsNull()) {
                decoder_template.CopyFrom(*decoder_projection, ontology.GetAllocator());
            }

            /* project decoder attributes from the ontology root */
            Value default_decoder(kObjectType);
            project_json_value(decoder_template, ontology, default_decoder, ontology);

            /* aseemble a default barcode */
            string barcode_projection_uri(key);
            barcode_projection_uri.append(":barcode");
            const Value* barcode_projection(find_projection(barcode_projection_uri));

            Value barcode_template(kObjectType);
            if(barcode_projection != NULL && !barcode_projection->IsNull()) {
                barcode_template.CopyFrom(*barcode_projection, ontology.GetAllocator());
            }

            /* project barcode attributes from the ontology root */
            Value default_barcode(kObjectType);
            project_json_value(barcode_template, ontology, default_barcode, ontology);

            int32_t index(0);
            if(reference->value.IsObject()) {
                try {
                    compile_decoder(reference->value, index, default_decoder, default_barcode);
                } catch(ConfigurationError& error) {
                    error.message.insert(0, string(key) + " decoder : ");
                    throw error;
                }
            } else if(reference->value.IsArray()) {
                try {
                    for(auto& element : reference->value.GetArray()) {
                        compile_decoder(element, index, default_decoder, default_barcode);
                    }
                } catch(ConfigurationError& error) {
                    error.message.insert(0, string(key) + " decoder at " + to_string(index) + " : ");
                    throw error;
                }
            }
            clean_json_value(reference->value, ontology);
        }
    }
};
void MultiplexJob::compile_decoder(Value& value, int32_t& index, const Value& default_decoder, const Value& default_barcode) {
    if(value.IsObject()) {
        encode_key_value("index", index, value, ontology);
        compile_codec(value, default_decoder, default_barcode);
        compile_decoder_transformation(value);
        ++index;
    }
};
void MultiplexJob::compile_codec(Value& value, const Value& default_decoder, const Value& default_barcode) {
    if(value.IsObject()) {
        /* overlay on top of the default decoder */
        merge_json_value(default_decoder, value, ontology);
        clean_json_value(value, ontology);

        /* compute default barcode induced by the codec */
        Value default_codec_barcode(kObjectType);
        project_json_value(default_barcode, value, default_codec_barcode, ontology);

        string buffer;
        int32_t barcode_index(0);
        double total_concentration(0);
        set< string > unique_barcode_id;

        double noise(decode_value_by_key< double >("noise", value));

        /* apply barcode default on undetermined barcode or create it from the default if one was not explicitly specified */
        Value::MemberIterator reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()){
            merge_json_value(default_codec_barcode, reference->value, ontology);
        } else {
            value.AddMember (
                Value("undetermined", ontology.GetAllocator()).Move(),
                Value(default_codec_barcode, ontology.GetAllocator()).Move(),
                ontology.GetAllocator()
            );
        }

        reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()){
            encode_key_value("index", barcode_index, reference->value, ontology);
            if(infer_ID("ID", buffer, reference->value, true)) {
                unique_barcode_id.emplace(buffer);
            }
            encode_key_value("concentration", noise, reference->value, ontology);
            ++barcode_index;
        }

        reference = value.FindMember("codec");
        if(reference != value.MemberEnd()){
            if(reference->value.IsObject()) {
                Value& codec(reference->value);

                for(auto& record : codec.GetObject()) {
                    merge_json_value(default_codec_barcode, record.value, ontology);
                    encode_key_value("index", barcode_index, record.value, ontology);
                    if(infer_ID("ID", buffer, record.value)) {
                        if(!unique_barcode_id.count(buffer)) {
                            unique_barcode_id.emplace(buffer);
                        } else {
                            string duplicate(record.name.GetString(), record.value.GetStringLength());
                            throw ConfigurationError("duplicate " + duplicate + " barcode");
                        }
                    }

                    double concentration(decode_value_by_key< double >("concentration", record.value));
                    if(concentration >= 0) {
                        total_concentration += concentration;
                    } else { throw ConfigurationError("barcode concentration must be a positive number");  }

                    ++barcode_index;
                }

                if(total_concentration > 0) {
                    const double factor((1.0 - noise) / total_concentration);
                    for(auto& record : codec.GetObject()) {
                        double concentration(decode_value_by_key< double >("concentration", record.value));
                        encode_key_value("concentration", concentration * factor, record.value, ontology);
                    }
                } else { throw ConfigurationError("total pool concentration is not a positive number"); }
            } else { throw ConfigurationError("codec element must be a dictionary"); }
        }
    }
};
void MultiplexJob::compile_decoder_transformation(Value& value) {
    if(value.HasMember("transform")) {
        compile_transformation(value);

        Rule rule(decode_value_by_key< Rule >("transform", value));
        int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));

        /* validate all tokens refer to an existing input segment */
        for(auto& token : rule.token_array) {
            if(!(token.input_segment_index < input_segment_cardinality)) {
                throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
            }
        }

        /* annotate the decoder with cardinality information from the transfortmation */
        int32_t nucleotide_cardinality(0);
        vector< int32_t > barcode_length(rule.output_segment_cardinality, 0);
        for(auto& transform : rule.transform_array) {
            if(transform.token.constant()) {
                if(!transform.token.empty()) {
                    barcode_length[transform.output_segment_index] += transform.token.length();
                    nucleotide_cardinality += transform.token.length();
                } else { throw ConfigurationError("multiplex barcode token " + string(transform.token) + " is empty"); }
            } else { throw ConfigurationError("barcode token " + string(transform.token) + " is not fixed width"); }
        }
        encode_key_value("segment cardinality", rule.output_segment_cardinality, value, ontology);
        encode_key_value("nucleotide cardinality", nucleotide_cardinality, value, ontology);
        encode_key_value("barcode length", barcode_length, value, ontology);

        /* annotate each barcode element with the barcode segment cardinality */
        Value::MemberIterator reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()) {
            if(!reference->value.IsNull()) {
                Value& undetermined(reference->value);

                /* explicitly define a null barcode segment for the right dimension in the undetermined */
                Value barcode(kArrayType);
                for(size_t i(0); i < barcode_length.size(); ++i) {
                    string sequence(barcode_length[i], '=');
                    barcode.PushBack(Value(sequence.c_str(), sequence.size(), ontology.GetAllocator()).Move(), ontology.GetAllocator());
                }
                undetermined.RemoveMember("barcode");
                undetermined.AddMember(Value("barcode", ontology.GetAllocator()).Move(), barcode.Move(), ontology.GetAllocator());
                encode_key_value("segment cardinality", rule.output_segment_cardinality, undetermined, ontology);
            }
        }

        reference = value.FindMember("codec");
        if(reference != value.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    for(auto& record : reference->value.GetObject()) {
                        encode_key_value("segment cardinality", rule.output_segment_cardinality, record.value, ontology);
                    }
                }
            }
        }

        CodecMetric metric(value);
        metric.compile_barcode_tolerance(value, ontology);

    }
};
void MultiplexJob::compile_output_transformation() {
    const int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));

    /* if the output transform was not defined add an empty dictionary */
    if(!ontology.HasMember("transform")) {
        ontology.AddMember("transform", Value(kObjectType).Move(), ontology.GetAllocator());
    }

    /* if transform does not define a token array route all input segments to the output verbatim */
    if(!ontology["transform"].HasMember("token")) {
        Value token_array(kArrayType);
        for(int32_t i(0); i < input_segment_cardinality; ++i) {
            string token(to_string(i) + "::");
            token_array.PushBack(Value(token.c_str(), token.size(), ontology.GetAllocator()), ontology.GetAllocator());
        }
        ontology["transform"].AddMember("token", token_array.Move(), ontology.GetAllocator());
    }
    compile_transformation(ontology);
};
void MultiplexJob::compile_output() {
    /* load output transform */
    compile_output_transformation();
    Rule rule(decode_value_by_key< Rule >("transform", ontology));

    const int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));
    for(auto& token : rule.token_array) {
        if(!(token.input_segment_index < input_segment_cardinality)) {
            throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
        }
    }
    const int32_t output_segment_cardinality(rule.output_segment_cardinality);
    encode_key_value("output segment cardinality", output_segment_cardinality, ontology, ontology);

    Platform platform(decode_value_by_key< Platform >("platform", ontology));
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    uint8_t phred_offset(decode_value_by_key< uint8_t >("output phred offset", ontology));

    Value::MemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            Value& value(reference->value);
            expand_url_value_by_key("base output url", value, ontology);
            URL base(decode_value_by_key< URL >("base output url", value));

            unordered_map< URL, unordered_map< int32_t, int > > feed_resolution;
            Value::MemberIterator reference = value.FindMember("undetermined");
            if(reference != value.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    int32_t index(decode_value_by_key< int32_t >("index", reference->value));
                    encode_key_value("TC", output_segment_cardinality, reference->value, ontology);
                    pad_url_array_by_key("output", reference->value, output_segment_cardinality);
                    expand_url_array_by_key("output", reference->value, ontology, IoDirection::OUT);
                    relocate_url_array_by_key("output", reference->value, ontology, base);

                    list< URL > feed_url_array;
                    if(decode_value_by_key< list< URL > >("output", feed_url_array, reference->value)) {
                        for(auto& url : feed_url_array) {
                            ++(feed_resolution[url][index]);
                        }
                    }
                }
            }
            reference = value.FindMember("codec");
            if(reference != value.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    if(reference->value.IsObject()) {
                        for(auto& record : reference->value.GetObject()) {
                            int32_t index(decode_value_by_key< int32_t >("index", record.value));
                            encode_key_value("TC", output_segment_cardinality, record.value, ontology);
                            pad_url_array_by_key("output", record.value, output_segment_cardinality);
                            expand_url_array_by_key("output", record.value, ontology, IoDirection::OUT);
                            relocate_url_array_by_key("output", record.value, ontology, base);

                            list< URL > feed_url_array;
                            if(decode_value_by_key< list< URL > >("output", feed_url_array, record.value)) {
                                for(auto& url : feed_url_array) {
                                    ++(feed_resolution[url][index]);
                                }
                            }
                        }
                    }
                }
            }

            if(feed_resolution.size() > 0) {
                unordered_map< URL, Value > feed_ontology_by_url;
                int32_t index(0);
                for(const auto& url_record : feed_resolution) {
                    const URL& url = url_record.first;
                    int resolution(0);
                    for(const auto& record : url_record.second) {
                        if(resolution == 0) {
                            resolution = record.second;
                        } else if(resolution != record.second) {
                            throw ConfigurationError("inconsistent resolution for " + string(url));
                        }
                    }
                    Value proxy(kObjectType);
                    encode_key_value("index", index, proxy, ontology);
                    encode_key_value("url", url, proxy, ontology);
                    encode_key_value("direction", IoDirection::OUT, proxy, ontology);
                    encode_key_value("platform", platform, proxy, ontology);
                    encode_key_value("capacity", buffer_capacity * resolution, proxy, ontology);
                    encode_key_value("resolution", resolution, proxy, ontology);
                    encode_key_value("phred offset", phred_offset, proxy, ontology);
                    feed_ontology_by_url.emplace(make_pair(url, move(proxy)));
                    ++index;
                }

                reference = value.FindMember("undetermined");
                if(reference != value.MemberEnd()) {
                    if(!reference->value.IsNull()) {
                        if(reference->value.IsObject()) {
                            list< URL > feed_url_array;
                            if(decode_value_by_key< list< URL > >("output", feed_url_array, reference->value)) {
                                Value feed_by_segment(kArrayType);
                                for(const auto& url : feed_url_array) {
                                    const Value& proxy(feed_ontology_by_url[url]);
                                    feed_by_segment.PushBack(Value(proxy, ontology.GetAllocator()).Move(), ontology.GetAllocator());
                                }
                                reference->value.RemoveMember("output feed by segment");
                                reference->value.AddMember("output feed by segment", feed_by_segment.Move(), ontology.GetAllocator());
                            }
                        }
                    }
                }

                reference = value.FindMember("codec");
                if(reference != value.MemberEnd()) {
                    if(!reference->value.IsNull()) {
                        if(reference->value.IsObject()) {
                            for(auto& record : reference->value.GetObject()) {
                                list< URL > feed_url_array;
                                if(decode_value_by_key< list< URL > >("output", feed_url_array, record.value)) {
                                    Value feed_by_segment(kArrayType);
                                    for(const auto& url : feed_url_array) {
                                        const Value& proxy(feed_ontology_by_url[url]);
                                        feed_by_segment.PushBack(Value(proxy, ontology.GetAllocator()).Move(), ontology.GetAllocator());
                                    }
                                    record.value.RemoveMember("output feed by segment");
                                    record.value.AddMember("output feed by segment", feed_by_segment.Move(), ontology.GetAllocator());
                                }
                            }
                        }
                    }
                }

                Value feed_array(kArrayType);
                for(auto& record : feed_ontology_by_url) {
                    feed_array.PushBack(record.second.Move(), ontology.GetAllocator());
                }
                ontology.RemoveMember("output feed");
                ontology.AddMember("output feed", feed_array.Move(), ontology.GetAllocator());
            }
            cross_validate_io();
        }
    }
};
void MultiplexJob::compile_transformation(Value& value) {
    /* add the default observation if one was not specificed.
       default observation will treat every token as a segment */
    if(value.IsObject()) {
        Value::MemberIterator reference = value.FindMember("transform");
        if(reference != value.MemberEnd()) {
            if(reference->value.IsObject()) {
                Value& transform_element(reference->value);
                reference = transform_element.FindMember("token");
                if(reference != transform_element.MemberEnd()) {
                    if(reference->value.IsArray()) {
                        int32_t token_cardinality(reference->value.Size());

                        reference = transform_element.FindMember("segment pattern");
                        if(reference == transform_element.MemberEnd() || reference->value.IsNull() || (reference->value.IsArray() && reference->value.Empty())) {
                            Value observation(kArrayType);
                            for(int32_t i(0); i < token_cardinality; ++i) {
                                string element(to_string(i));
                                observation.PushBack(Value(element.c_str(), element.size(),ontology.GetAllocator()).Move(), ontology.GetAllocator());
                            }
                            transform_element.RemoveMember("segment pattern");
                            transform_element.AddMember(Value("segment pattern", ontology.GetAllocator()).Move(), observation.Move(), ontology.GetAllocator());
                        }
                    } else { throw ConfigurationError("transform token element is not an array"); }
                } else { throw ConfigurationError("transform element is missing a token array"); }
            }
        }
    }
};
bool MultiplexJob::infer_PU(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
    buffer.clear();
    string suffix;
    if(!decode_value_by_key< string >(key, suffix, container)) {
        if(!undetermined) {
            list< string > barcode;
            if(decode_value_by_key< list< string > >("barcode", barcode, container)) {
                for(auto& segment : barcode) {
                    suffix.append(segment);
                }
            }
        } else { suffix.assign("undetermined"); }

        if(!suffix.empty()) {
            if(decode_value_by_key< string >("flowcell id", buffer, container)) {
                buffer.push_back(':');
                int32_t lane;
                if(decode_value_by_key< int32_t >("flowcell lane number", lane, container)) {
                    buffer.append(to_string(lane));
                    buffer.push_back(':');
                }
            }
            buffer.append(suffix);
            encode_key_value(key, buffer, container, ontology);
            return true;
        } else { return false; }
    } else { return true; }
};
bool MultiplexJob::infer_ID(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
    buffer.clear();
    if(!decode_value_by_key< string >(key, buffer, container)) {
        if(infer_PU("PU", buffer, container, undetermined)) {
            encode_key_value(key, buffer, container, ontology);
            return true;
        } else { return false; }
    } else { return true; }
};
void MultiplexJob::pad_url_array_by_key(const Value::Ch* key, Value& container, const int32_t& cardinality) {
    list< URL > array;
    if(decode_value_by_key< list< URL > >(key, array, container)) {
        if(!array.empty()) {
            if(static_cast< int32_t >(array.size()) != cardinality) {
                if(array.size() == 1) {
                    while(static_cast< int32_t >(array.size()) < cardinality) {
                        array.push_back(array.front());
                    }
                    encode_key_value(key, array, container, ontology);
                } else { throw ConfigurationError("incorrect number of output URLs in channel"); }
            }
        }
    }
};
void MultiplexJob::cross_validate_io() {
    list< URL > input_array(decode_value_by_key< list< URL > >("input", ontology));

    set< URL > input;
    for(auto& url : input_array) {
        input.emplace(url);
    }

    Value::MemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            reference = reference->value.FindMember("codec");
            if(reference != ontology.MemberEnd()) {
                if(reference->value.IsObject()) {
                    for(auto& record : reference->value.GetObject()) {
                        list< URL > output;
                        if(decode_value_by_key< list< URL > >("output", output, record.value)) {
                            for(auto& url : output) {
                                if(input.count(url) > 0) {
                                    throw ConfigurationError("URL " + string(url) + " is used for both input and output");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};
void MultiplexJob::validate_decoder_group(const Value::Ch* key) {
    Value::MemberIterator reference = ontology.FindMember(key);
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                validate_decoder(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& decoder : reference->value.GetArray()) {
                    if(!decoder.IsNull()) {
                        validate_decoder(decoder);
                    }
                }
            }
        }
    }
};
void MultiplexJob::validate_decoder(Value& value) {
    if(value.IsObject()) {
        if(value.HasMember("codec")) {
            double confidence_threshold;
            if(decode_value_by_key< double >("confidence threshold", confidence_threshold, value)) {
                if(confidence_threshold < 0 || confidence_threshold > 1) {
                    throw ConfigurationError("confidence threshold value " + to_string(confidence_threshold) + " not between 0 and 1");
                }
            }

            double noise;
            if(decode_value_by_key< double >("noise", noise, value)) {
                if(noise < 0 || noise > 1) {
                    throw ConfigurationError("noise value " + to_string(noise) + " not between 0 and 1");
                }
            }
        }
    }
};
void MultiplexJob::validate_url_accessibility() {
    URL url;
    Value::MemberIterator reference = ontology.FindMember("input feed");
    if(reference != ontology.MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_readable()) {
                    throw IOError("could not open " + string(url) + " for reading");
                }
            }
        }
    }

    reference = ontology.FindMember("output feed");
    if(reference != ontology.MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_writable()) {
                    throw IOError("could not open " + string(url) + " for writing");
                }
            }
        }
    }
};
void MultiplexJob::load_thread_pool() {
    int32_t threads(decode_value_by_key< int32_t >("threads", ontology));
    thread_pool.pool = hts_tpool_init(threads);
    if(!thread_pool.pool) { throw InternalError("error creating thread pool"); }
};
void MultiplexJob::load_input() {
    if(input_feed_by_index.empty()) {
        /*  Decode feed_proxy_array, a local list of input feed proxy.
            The list has already been enumerated by the interface
            and contains only unique url references
        */
        list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("input feed", ontology));

        /*  Initialized the hfile reference and verify input format */
        for(auto& proxy : feed_proxy_array) {
            proxy.probe();
        };

        /*  Load feed_by_url, a local map of input feeds by url, from the proxy.
            Populate input_feed_by_index used to enumerate threaded access to the input feeds */
        unordered_map< URL, Feed* > feed_by_url(feed_proxy_array.size());
        for(auto& proxy : feed_proxy_array) {
            Feed* feed(NULL);
            switch(proxy.kind()) {
                case FormatKind::FASTQ: {
                    feed = new FastqFeed(proxy);
                    break;
                };
                case FormatKind::HTS: {
                    feed = new HtsFeed(proxy);
                    break;
                };
                case FormatKind::DEV_NULL: {
                    feed = new NullFeed(proxy);
                    break;
                };
                default: {
                    throw InternalError("unknown input format " + string(proxy.url));
                    break;
                };
            }
            feed->set_thread_pool(&thread_pool);
            input_feed_by_index.push_back(feed);
            feed_by_url.emplace(make_pair(proxy.url, feed));
        }

        list< URL > url_by_segment(decode_value_by_key< list < URL > >("input", ontology));

        /*  Populate the input_feed_by_segment array */
        input_feed_by_segment.reserve(url_by_segment.size());
        for(auto& url : url_by_segment) {
            const auto& record = feed_by_url.find(url);
            if(record != feed_by_url.end()) {
                input_feed_by_segment.push_back(record->second);
            } else {
                throw InternalError("missing feed for URL " + string(url) + " referenced in input proxy segment array");
            }
        }
    }
};
void MultiplexJob::load_output() {
    /*  Decode feed_proxy_array, a local list of output feed proxy.
        The list has already been enumerated by the environment
        and contains only unique url references
    */
    list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("output feed", ontology));
    HeadPGAtom program(decode_value_by_key< HeadPGAtom >("program", ontology));

    /*  Register the read group elements on the feed proxy so it can be added to SAM header
        if a URL is present in the channel output that means the channel writes output to that file
        and the read group should be added to the header of that file.
    */
    map< URL, FeedProxy* > feed_proxy_by_url;
    for(auto& proxy : feed_proxy_array) {
        proxy.register_pg(program);
        feed_proxy_by_url.emplace(make_pair(proxy.url, &proxy));
    };

    Value::ConstMemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                const Value& multiplex(reference->value);

                reference = multiplex.FindMember("undetermined");
                if(reference != multiplex.MemberEnd()) {
                    HeadRGAtom rg(reference->value);
                    list< URL > output(decode_value_by_key< list< URL > >("output", reference->value));
                    for(auto& url : output) {
                        feed_proxy_by_url[url]->register_rg(rg);
                    }
                }
                reference = multiplex.FindMember("codec");
                if(reference != multiplex.MemberEnd()) {
                    for(auto& record : reference->value.GetObject()) {
                        HeadRGAtom rg(record.value);
                        list< URL > output(decode_value_by_key< list< URL > >("output", record.value));
                        for(auto& url : output) {
                            feed_proxy_by_url[url]->register_rg(rg);
                        }
                    }
                }
            } else { throw ConfigurationError("multiplex element must be a dictionary"); }
        }
    }

    /*  Initialized the hfile reference */
    for(auto& proxy : feed_proxy_array) {
        proxy.probe();
    };

    /*  Load feed_by_url, a local map of output feeds by url, from the proxy.
        Populate output_feed_by_index used to enumerate threaded access to the output feeds */
    output_feed_by_url.reserve(feed_proxy_array.size());
    for(auto& proxy : feed_proxy_array) {
            Feed* feed(NULL);
            switch(proxy.kind()) {
                case FormatKind::FASTQ: {
                    feed = new FastqFeed(proxy);
                    break;
                };
                case FormatKind::HTS: {
                    feed = new HtsFeed(proxy);
                    break;
                };
                case FormatKind::DEV_NULL: {
                    feed = new NullFeed(proxy);
                    break;
                };
                default: {
                    throw InternalError("unknown output format " + string(proxy.url));
                    break;
                };
            }
            feed->set_thread_pool(&thread_pool);
            output_feed_by_index.push_back(feed);
            output_feed_by_url.emplace(make_pair(proxy.url, feed));
    }
};
void MultiplexJob::load_pivot() {
    int32_t threads(decode_value_by_key< int32_t >("threads", ontology));
    for(int32_t index(0); index < threads; ++index) {
        pivot_array.emplace_back(*this, index);
    }
};
void MultiplexJob::populate_channel(Channel& channel) {
    map< int32_t, Feed* > feed_by_index;
    channel.output_feed_by_segment.reserve(channel.output_feed_url_by_segment.size());
    for(const auto& url : channel.output_feed_url_by_segment) {
        Feed* feed(output_feed_by_url[url]);
        channel.output_feed_by_segment.emplace_back(feed);
        if(feed_by_index.count(feed->index) == 0) {
            feed_by_index.emplace(make_pair(feed->index, feed));
        }
    }
    channel.output_feed_by_segment.shrink_to_fit();

    channel.output_feed_lock_order.reserve(feed_by_index.size());
    for(auto& record : feed_by_index) {
        /* /dev/null is not really being written to so we don't need to lock it */
        if(!record.second->is_dev_null()) {
            channel.output_feed_lock_order.push_back(record.second);
        }
    }
    channel.output_feed_lock_order.shrink_to_fit();
};
void MultiplexJob::print_global_instruction(ostream& o) const {
    o << setprecision(16);
    o << "Environment " << endl << endl;
    // o << "    Version                                     " << interface.application_version << endl;

    URL base_input_url;
    decode_value_by_key< URL >("base input url", base_input_url, ontology);
    o << "    Base input URL                              " << base_input_url << endl;

    URL base_output_url;
    decode_value_by_key< URL >("base input url", base_output_url, ontology);
    o << "    Base output URL                             " << base_output_url << endl;

    Platform platform;
    decode_value_by_key< Platform >("platform", platform, ontology);
    o << "    Platform                                    " << platform << endl;

    bool disable_quality_control;
    decode_value_by_key< bool >("disable quality control", disable_quality_control, ontology);
    o << "    Quality tracking                            " << (disable_quality_control ? "disabled" : "enabled") << endl;

    bool include_filtered;
    decode_value_by_key< bool >("include filtered", include_filtered, ontology);
    o << "    Include non PF reads                        " << (include_filtered ? "enabled" : "disabled") << endl;

    uint8_t input_phred_offset;
    decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, ontology);
    o << "    Input Phred offset                          " << to_string(input_phred_offset) << endl;

    uint8_t output_phred_offset;
    decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, ontology);
    o << "    Output Phred offset                         " << to_string(output_phred_offset) << endl;

    int32_t leading_segment_index;
    decode_value_by_key< int32_t >("leading segment index", leading_segment_index, ontology);
    o << "    Leading segment index                       " << to_string(leading_segment_index) << endl;

    int32_t buffer_capacity;
    decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, ontology);
    o << "    Feed buffer capacity                        " << to_string(buffer_capacity) << endl;

    int32_t threads;
    decode_value_by_key< int32_t >("threads", threads, ontology);
    o << "    Threads                                     " << to_string(threads) << endl;
    o << endl;
};
void MultiplexJob::print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const {
    Value::ConstMemberIterator reference = ontology.FindMember(key);
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            o << head << endl << endl;
            if(reference->value.IsObject()) {
                print_codec_instruction(reference->value, false, o);
            } else if(reference->value.IsArray()) {
                bool plural(reference->value.Size() > 1);
                for(auto& decoder : reference->value.GetArray()) {
                    if(!decoder.IsNull()) {
                        print_codec_instruction(decoder, plural, o);
                    }
                }
            }
        }
    }
};
void MultiplexJob::print_codec_instruction(const Value& value, const bool& plural, ostream& o) const {
    if(!value.IsNull()) {
        if(plural) {
            int32_t index;
            if(decode_value_by_key< int32_t >("index", index, value)) {
                o << "  Decoder No." << to_string(index) << endl << endl;
            }
        }

        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
        o << "    Decoding algorithm                          " << algorithm << endl;

        uint8_t quality_masking_threshold;
        if(decode_value_by_key< uint8_t >("quality masking threshold", quality_masking_threshold, value)) {
            if(quality_masking_threshold > 0) {
                o << "    Quality masking threshold                   " << int32_t(quality_masking_threshold) << endl;
            }
        }

        vector< int32_t > shannon_bound;
        if(decode_value_by_key< vector< int32_t > >("shannon bound", shannon_bound, value)) {
            o << "    Shannon bound                              ";
            for(auto& element : shannon_bound) { o << " " << element; }
            o << endl;
        }

        if(algorithm == Algorithm::MDD) {
            vector< int32_t > distance_tolerance;
            if(decode_value_by_key< vector< int32_t > >("distance tolerance", distance_tolerance, value)) {
                o << "    Distance tolerance                          ";
                for(auto& element : distance_tolerance) { o << " " << element; }
                o << endl;
            }
        }

        if(algorithm == Algorithm::PAMLD) {
            double noise(decode_value_by_key< double >("noise", value));
            o << "    Noise                                       " << noise << endl;

            double confidence_threshold(decode_value_by_key< double >("confidence threshold", value));
            o << "    Confidence threshold                        " << confidence_threshold << endl;
        }

        int32_t segment_cardinality(decode_value_by_key< int32_t >("segment cardinality", value));
        if(segment_cardinality > 0) {
            o << "    Segment cardinality                         " << to_string(segment_cardinality) << endl;

            int32_t nucleotide_cardinality;
            if(decode_value_by_key< int32_t >("nucleotide cardinality", nucleotide_cardinality, value)) {
                o << "    Nucleotide cardinality                      " << to_string(nucleotide_cardinality) << endl;
            }

            if(segment_cardinality > 1) {
                vector< int32_t > barcode_length;
                if(decode_value_by_key< vector< int32_t > >("barcode length", barcode_length, value)) {
                    o << "    Barcode segment length                      ";
                    for(const auto& v : barcode_length) {
                        o << to_string(v) << " ";
                    }
                    o << endl;
                }
            }

            o << endl << "    Transform" << endl;

            if(ontology.HasMember("transform")) {
                Rule rule(decode_value_by_key< Rule >("transform", value));
                o << endl;
                for(auto& token : rule.token_array) {
                    o << "        Token No." << token.index << endl;
                    o << "            Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
                    o << "            Pattern       " << string(token) << endl;
                    o << "            Description   ";
                    o << token.description() << endl;
                    o << endl;
                }
                o << "        Assembly instruction" << endl;
                for(const auto& transform : rule.transform_array) {
                    o << "            " << transform.description() << endl;
                }
                o << endl;

                if(display_distance()) {
                    CodecMetric metric(value);
                    metric.describe(o);
                }
            }
        }
        o << endl;

        Value::ConstMemberIterator reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()) {
            const string key(reference->name.GetString(), reference->name.GetStringLength());
            print_channel_instruction(key, reference->value, o);
        }
        reference = value.FindMember("codec");
        if(reference != value.MemberEnd()) {
            if(reference->value.IsObject()) {
                for(const auto& record : reference->value.GetObject()) {
                    const string key(record.name.GetString(), record.name.GetStringLength());
                    print_channel_instruction(key, record.value, o);
                }
            }
        }
    }
};
void MultiplexJob::print_channel_instruction(const string& key, const Value& value, ostream& o) const {
    if(value.IsObject()) {
        o << "    Barcode " << key << endl;

        string buffer;
        if(decode_value_by_key< string >("ID", buffer, value)) { o << "        ID : " << buffer << endl; }
        if(decode_value_by_key< string >("PU", buffer, value)) { o << "        PU : " << buffer << endl; }
        if(decode_value_by_key< string >("LB", buffer, value)) { o << "        LB : " << buffer << endl; }
        if(decode_value_by_key< string >("SM", buffer, value)) { o << "        SM : " << buffer << endl; }
        if(decode_value_by_key< string >("DS", buffer, value)) { o << "        DS : " << buffer << endl; }
        if(decode_value_by_key< string >("DT", buffer, value)) { o << "        DT : " << buffer << endl; }
        if(decode_value_by_key< string >("PL", buffer, value)) { o << "        PL : " << buffer << endl; }
        if(decode_value_by_key< string >("PM", buffer, value)) { o << "        PM : " << buffer << endl; }
        if(decode_value_by_key< string >("CN", buffer, value)) { o << "        CN : " << buffer << endl; }
        if(decode_value_by_key< string >("FO", buffer, value)) { o << "        FO : " << buffer << endl; }
        if(decode_value_by_key< string >("KS", buffer, value)) { o << "        KS : " << buffer << endl; }
        if(decode_value_by_key< string >("PI", buffer, value)) { o << "        PI : " << buffer << endl; }
        if(decode_value_by_key< string >("FS", buffer, value)) { o << "        FS : " << buffer << endl; }
        if(decode_value_by_key< string >("CO", buffer, value)) { o << "        CO : " << buffer << endl; }

        int32_t index(decode_value_by_key< int32_t >("index", value));
        if(index > 0) {
            double concentration;
            if(decode_value_by_key< double >("concentration", concentration, value)) {
                o << "        Concentration : " << concentration << endl;
            }

            Barcode barcode(value);
            if(!barcode.empty()) { o << "        Barcode       : " << barcode.iupac_ambiguity() << endl; }
        }

        int32_t segment_index(0);
        list< URL > output;
        if(decode_value_by_key< list< URL > >("output", output, value)) {
            for(auto& url : output) {
                o << "        Segment No." + to_string(segment_index) + "  : " << url << endl;
                ++segment_index;
            }
        }
        o << endl;
    }
};
void MultiplexJob::print_feed_instruction(const Value::Ch* key, ostream& o) const {
    Value::ConstMemberIterator reference = ontology.FindMember(key);
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                if(!reference->value.Empty()) {
                    for(const auto& element : reference->value.GetArray()) {
                        IoDirection direction(decode_value_by_key< IoDirection >("direction", element));
                        int32_t index(decode_value_by_key< int32_t >("index", element));
                        int32_t resolution(decode_value_by_key< int32_t >("resolution", element));
                        int32_t capacity(decode_value_by_key< int32_t >("capacity", element));
                        URL url(decode_value_by_key< URL >("url", element));
                        Platform platform(decode_value_by_key< Platform >("platform", element));
                        uint8_t phred_offset(decode_value_by_key< uint8_t >("phred offset", element));

                        o << "    ";
                        switch (direction) {
                            case IoDirection::IN:
                                o << "Input";
                                break;
                            case IoDirection::OUT:
                                o << "Output";
                                break;
                            default:
                                break;
                        }
                        o << " feed No." << index << endl;
                        o << "        Type : " << url.type() << endl;
                        o << "        Resolution : " << resolution << endl;
                        o << "        Phred offset : " << to_string(phred_offset) << endl;
                        o << "        Platform : " << platform << endl;
                        o << "        Buffer capacity : " << capacity << endl;
                        o << "        URL : " << url << endl;
                        o << endl;
                    }
                }
            }
        }
    }
};
void MultiplexJob::print_input_instruction(ostream& o) const {
    o << "Input " << endl << endl;

    int32_t input_segment_cardinality;
    if(decode_value_by_key< int32_t >("input segment cardinality", input_segment_cardinality, ontology)) {
        o << "    Input segment cardinality                   " << to_string(input_segment_cardinality) << endl;
    }

    list< URL > input_url_array;
    if(decode_value_by_key< list< URL > >("input", input_url_array, ontology)) {
        o << endl;
        int32_t url_index(0);
        for(auto& url : input_url_array) {
            o << "    Input segment No." << url_index << " : " << url << endl;
            ++url_index;
        }
        o << endl;
    }
    print_feed_instruction("input feed", o);
};
void MultiplexJob::print_transform_instruction(ostream& o) const {
    o << "Output transform" << endl << endl;

    int32_t output_segment_cardinality;
    if(decode_value_by_key< int32_t >("output segment cardinality", output_segment_cardinality, ontology)) {
        o << "    Output segment cardinality                  " << to_string(output_segment_cardinality) << endl;
    }

    Rule rule(decode_value_by_key< Rule >("transform", ontology));
    o << endl;
    for(auto& token : rule.token_array) {
        o << "    Token No." << token.index << endl;
        o << "        Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
        o << "        Pattern       " << string(token) << endl;
        o << "        Description   ";
        o << token.description() << endl;
        o << endl;
    }
    o << "    Assembly instruction" << endl;
    for(const auto& transform : rule.transform_array) {
        o << "        " << transform.description() << endl;
    }
    o << endl;
};
void MultiplexJob::print_multiplex_instruction(ostream& o) const {
    print_codec_group_instruction("multiplex", "Mutliplex decoding", o);
    print_feed_instruction("output feed", o);
};
void MultiplexJob::print_molecular_instruction(ostream& o) const {
    print_codec_group_instruction("molecular", "Molecular decoding", o);
};
void MultiplexJob::print_cellular_instruction(ostream& o) const {
    print_codec_group_instruction("cellular", "Cellular decoding", o);
};

MultiplexPivot::MultiplexPivot(MultiplexJob& job, const int32_t& index) try :
    index(index),
    platform(decode_value_by_key< Platform >("platform", job.ontology)),
    leading_segment_index(decode_value_by_key< int32_t >("leading segment index", job.ontology)),
    input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", job.ontology)),
    output_segment_cardinality(decode_value_by_key< int32_t >("output segment cardinality", job.ontology)),
    input(input_segment_cardinality, platform, leading_segment_index),
    output(output_segment_cardinality, platform, leading_segment_index),
    multiplex(NULL),
    input_accumulator(job.ontology),
    output_accumulator(find_value_by_key("multiplex", job.ontology)),
    job(job),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", job.ontology)),
    template_rule(decode_value_by_key< Rule >("transform", job.ontology)) {

    load_multiplex_decoding();
    load_molecular_decoding();
    load_cellular_decoding();
    clear();

    } catch(Error& error) {
        error.push("MultiplexPivot");
        throw;
};
void MultiplexPivot::load_multiplex_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("multiplex");
    if(reference != job.ontology.MemberEnd()) {
        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", reference->value));
        switch (algorithm) {
            case Algorithm::PAMLD: {
                MultiplexPAMLDecoder* pamld_decoder(new MultiplexPAMLDecoder(reference->value));
                pamld_decoder->unclassified.populate(job.output_feed_by_url);
                for(auto& channel : pamld_decoder->element_by_index) {
                    channel.populate(job.output_feed_by_url);
                }
                multiplex = pamld_decoder;
                break;
            };
            case Algorithm::MDD: {
                MultiplexMDDecoder* mdd_decoder(new MultiplexMDDecoder(reference->value));
                mdd_decoder->unclassified.populate(job.output_feed_by_url);
                for(auto& channel : mdd_decoder->element_by_index) {
                    channel.populate(job.output_feed_by_url);
                }
                multiplex = mdd_decoder;
                break;
            };
            case Algorithm::TRANSPARENT: {
                TransparentDecoder< Channel >* transparent_decoder(new TransparentDecoder< Channel >(reference->value));
                transparent_decoder->unclassified.populate(job.output_feed_by_url);
                multiplex = transparent_decoder;
                break;
            };
            default:
                throw ConfigurationError("unknown multiplex decoder algorithm");
                break;
        }
    }
};
void MultiplexPivot::load_molecular_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("molecular");
    if(reference != job.ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            molecular.reserve(1);
            load_molecular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            molecular.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_molecular_decoder(element);
            }
        }
    }
    molecular.shrink_to_fit();
};
void MultiplexPivot::load_molecular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch (algorithm) {
        case Algorithm::NAIVE: {
            MolecularNaiveDecoder* naive_decoder(new MolecularNaiveDecoder(value));
            molecular.emplace_back(naive_decoder);
            break;
        };
        default:
            break;
    }
};
void MultiplexPivot::load_cellular_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("cellular");
    if(reference != job.ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            molecular.reserve(1);
            load_cellular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            molecular.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_cellular_decoder(element);
            }
        }
    }
    cellular.shrink_to_fit();
};
void MultiplexPivot::load_cellular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch (algorithm) {
        case Algorithm::PAMLD: {
            CellularPAMLDecoder* paml_decoder(new CellularPAMLDecoder(value));
            cellular.emplace_back(paml_decoder);
            break;
        };
        case Algorithm::MDD: {
            CellularMDDecoder* md_decoder(new CellularMDDecoder(value));
            cellular.emplace_back(md_decoder);
            break;
        };
        default:
            break;
    }
};
