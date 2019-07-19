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

/*
    filter_incoming_qc_fail:    should we process incoming qc_fail reads?
    filter_outgoing_qc_fail:    should we write qc_fail reads?
    inherit_qc_fail:            should outgoing reads inherit qc_fail flag from incoming?
*/
#include "transcode.h"

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

Transcode::Transcode(Document& operation) try :
    Job(operation),
    count(0),
    pf_count(0),
    pf_fraction(0),
    end_of_input(false),
    decoded_nucleotide_cardinality(0),
    thread_pool({NULL, 0}) {

    } catch(Error& error) {
        error.push("Transcode");
        throw;
};
Transcode::~Transcode() {
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
bool Transcode::pull(Read& read) {
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

    /* update input counters */
    if(!end_of_input) {
        ++count;
        if(!read.qcfail()) {
            ++pf_count;
        }
    }

    /* release the locks on the input feeds in reverse order */
    for(auto feed_lock(feed_locks.rbegin()); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }
    return !end_of_input;
};
Transcode& Transcode::operator+=(const TranscodePivot& pivot) {
    if(sample_classifier != NULL) {
        *sample_classifier += *(pivot.sample_classifier);
    }
    if(!molecular_classifier_array.empty()) {
        for(size_t index(0); index < molecular_classifier_array.size(); ++index) {
            *(molecular_classifier_array[index]) += *(pivot.molecular_classifier_array[index]);
        }
    }
    if(!cellular_classifier_array.empty()) {
        for(size_t index(0); index < cellular_classifier_array.size(); ++index) {
            *(cellular_classifier_array[index]) += *(pivot.cellular_classifier_array[index]);
        }
    }
    return *this;
};

/* assemble */
void Transcode::assemble() {
    Job::assemble();
    apply_inheritance();
    clean_json_object(instruction, instruction);
};
void Transcode::apply_inheritance() {
    /* apply inheritence between decoders defined in the repository */
    apply_repository_inheritence("decoder", instruction, instruction);

    /* Apply inheritence in the indivitual classification categories */
    apply_topic_inheritance("multiplex");
    apply_topic_inheritance("molecular");
    apply_topic_inheritance("cellular");

    if(instruction.HasMember("transform")) {
        if(!instruction.HasMember("template")) {
            instruction.AddMember("template", Value(kObjectType).Move(), instruction.GetAllocator());
        }
        if(!instruction["template"].HasMember("transform")) {
            instruction["template"].AddMember("transform", Value(kObjectType).Move(), instruction.GetAllocator());
        }
        merge_json_value(instruction["transform"], instruction["template"]["transform"], instruction);
    }

    /* Remove the decoder repository, it is no longer needed in the compiled instruction */
    instruction.RemoveMember("decoder");

    sort_json_value(instruction, instruction);
};
void Transcode::apply_repository_inheritence(const Value::Ch* key, Value& container, Document& document) {
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
void Transcode::apply_topic_inheritance(const Value::Ch* key) {
    if(instruction.IsObject()) {
        Value::MemberIterator reference = instruction.FindMember(key);
        if(reference != instruction.MemberEnd()) {
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
void Transcode::apply_decoder_inheritance(Value& value) {
    if(value.IsObject()) {
        string base;
        Value::ConstMemberIterator reference;
        const Pointer decoder_repository_query("/decoder");
        if(decode_value_by_key< string >("base", base, value)) {
            const Value* decoder_repository(decoder_repository_query.Get(instruction));
            if(decoder_repository != NULL) {
                reference = decoder_repository->FindMember(base.c_str());
                if(reference != decoder_repository->MemberEnd()) {
                    merge_json_value(reference->value, value, instruction);
                } else { throw ConfigurationError("reference to an unknown base " + base); }
            }
        }
        value.RemoveMember("base");
        clean_json_value(value, instruction);
    }
};

/* compile */
void Transcode::compile() {
    ontology.CopyFrom(instruction, ontology.GetAllocator());
    remove_disabled_from_json_value(ontology, ontology);
    clean_json_object(ontology, ontology);

    ontology.RemoveMember("feed");
    ontology.RemoveMember("input segment cardinality");
    ontology.RemoveMember("output segment cardinality");
    ontology.RemoveMember("program");

    /* overlay on top of the default configuration */
    apply_default_ontology(ontology);

    /* overlay interactive parameters on top of the configuration */
    apply_interactive_ontology(ontology);

    /* compile a PG SAM header ontology with details about pheniqs */
    compile_PG();

    ontology.AddMember("feed", Value(kObjectType).Move(), ontology.GetAllocator());
    compile_input();
    compile_barcode_decoding();
    compile_output();

    compile_thread_model();
    clean_json_object(ontology, ontology);
    validate();
};
void Transcode::compile_PG() {
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
void Transcode::compile_input() {
    int32_t total_threads(decode_value_by_key< int32_t >("threads", ontology));

    int32_t htslib_threads;
    if(!decode_value_by_key("htslib threads", htslib_threads, ontology)) {
        htslib_threads = max(1, total_threads);
        encode_key_value("htslib threads", htslib_threads, ontology, ontology);
    }

    /* Populate the input_feed_by_index and input_feed_by_segment arrays */
    standardize_url_value_by_key("base input url", ontology, ontology, IoDirection::IN);
    URL base(decode_value_by_key< URL >("base input url", ontology));

    standardize_url_array_by_key("input", ontology, ontology, IoDirection::IN);
    relocate_url_array_by_key("input", ontology, ontology, base);

    /* Collect query parameters from all reoccurring references to the same path */
    unordered_map< string, URL > url_by_path;
    list< URL > feed_url_by_index(decode_value_by_key< list< URL > >("input", ontology));
    for(auto& url : feed_url_by_index) {
        auto record = url_by_path.find(url.path());
        if(record == url_by_path.end()) {
            url_by_path[url.path()] = url;
        } else {
            record->second.override_query(url);
        }
    }
    for(auto& url : feed_url_by_index) {
        url = url_by_path[url.path()];
    }
    encode_key_value("input", feed_url_by_index, ontology, ontology);

    if(is_sense_input_layout()) {
        compile_sensed_input();
    } else {
        compile_explicit_input();
    }

    /*  validate leading_segment_index */
    int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));
    int32_t leading_segment_index(decode_value_by_key< int32_t >("leading segment index", ontology));
    if(leading_segment_index >= input_segment_cardinality) {
        throw ConfigurationError("leading segment index " + to_string(leading_segment_index) + " references non existing input segment");
    }
};
void Transcode::compile_sensed_input() {
    int32_t feed_index(0);
    int32_t input_segment_cardinality(0);
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    uint8_t input_phred_offset(decode_value_by_key< uint8_t >("input phred offset", ontology));
    list< URL > explicit_url_array(decode_value_by_key< list< URL > >("input", ontology));
    Platform platform(decode_value_by_key< Platform >("platform", ontology));

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

    load_thread_pool();
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
        proxy.open();
        Segment segment;
        string feed_read_id;
        switch(proxy.kind()) {
            case FormatKind::FASTQ: {
                FastqFeed* fastq_feed = new FastqFeed(proxy);
                fastq_feed->set_thread_pool(&thread_pool);
                fastq_feed->initiate();
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
                hts_feed->initiate();
                hts_feed->replenish();
                if(hts_feed->peek(segment, resolution)) {
                    feed_read_id.assign(segment.name.s, segment.name.l);
                    resolution = segment.total_segments();
                }
                hts_feed->calibrate_resolution(resolution);
                /*
                    const HtsHead& head = ((HtsFeed*)feed)->get_header();
                    for(const auto& record : head.read_group_by_id) {}
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
    encode_key_value("input feed", input_feed_by_index, ontology["feed"], ontology);
    encode_key_value("input feed by segment", input_feed_by_segment, ontology["feed"], ontology);
};
void Transcode::compile_explicit_input() {
    int32_t feed_index(0);
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    uint8_t input_phred_offset(decode_value_by_key< uint8_t >("input phred offset", ontology));
    list< URL > explicit_url_array(decode_value_by_key< list< URL > >("input", ontology));
    Platform platform(decode_value_by_key< Platform >("platform", ontology));

    /* encode the input segment cardinality */
    encode_key_value("input segment cardinality", static_cast< int32_t>(explicit_url_array.size()), ontology, ontology);

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
    ontology["feed"].RemoveMember("input feed by segment");
    ontology["feed"].AddMember("input feed by segment", feed_by_segment_array.Move(), ontology.GetAllocator());

    Value feed_by_index_array(kArrayType);
    for(const auto& url : feed_url_by_index) {
        feed_by_index_array.PushBack(feed_ontology_by_url[url].Move(), ontology.GetAllocator());
    }
    feed_ontology_by_url.clear();

    ontology["feed"].RemoveMember("input feed");
    ontology["feed"].AddMember("input feed", feed_by_index_array.Move(), ontology.GetAllocator());
};
void Transcode::compile_transformation(Value& value) {
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
                        size_t token_cardinality(reference->value.Size());

                        reference = transform_element.FindMember("segment pattern");
                        if(reference == transform_element.MemberEnd() || reference->value.IsNull() || (reference->value.IsArray() && reference->value.Empty())) {
                            Value observation(kArrayType);
                            for(size_t i(0); i < token_cardinality; ++i) {
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
void Transcode::compile_barcode_decoding() {
    compile_topic("multiplex");
    compile_topic("molecular");
    compile_topic("cellular");
};
void Transcode::compile_topic(const Value::Ch* key) {
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
void Transcode::compile_decoder(Value& value, int32_t& index, const Value& default_decoder, const Value& default_barcode) {
    if(value.IsObject()) {
        encode_key_value("index", index, value, ontology);
        compile_codec(value, default_decoder, default_barcode);
        ++index;
    }
};
void Transcode::compile_codec(Value& value, const Value& default_decoder, const Value& default_barcode) {
    if(value.IsObject()) {
        /* overlay on top of the default decoder */
        merge_json_value(default_decoder, value, ontology);
        clean_json_value(value, ontology);

        /* compute default barcode induced by the codec */
        Value default_codec_barcode(kObjectType);
        project_json_value(default_barcode, value, default_codec_barcode, ontology);

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

        compile_decoder_transformation(value);

        string buffer;
        int32_t barcode_index(0);
        double total_concentration(0);
        set< string > unique_barcode_id;
        double noise(decode_value_by_key< double >("noise", value));

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
                            string duplicate(record.name.GetString(), record.name.GetStringLength());
                            throw ConfigurationError("duplicate " + duplicate + " barcode");
                        }
                    }
                    double concentration(decode_value_by_key< double >("concentration", record.value));
                    if(concentration >= 0) {
                        total_concentration += concentration;
                    } else { throw ConfigurationError("barcode concentration must be a positive number");  }
                    ++barcode_index;
                }

                int32_t nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", value));
                encode_key_value("barcode cardinality", barcode_index, value, ontology);
                decoded_nucleotide_cardinality += barcode_index * nucleotide_cardinality;

                if(total_concentration > 0) {
                    const double factor((1.0 - noise) / total_concentration);
                    for(auto& record : codec.GetObject()) {
                        double concentration(decode_value_by_key< double >("concentration", record.value));
                        encode_key_value("concentration", concentration * factor, record.value, ontology);
                    }
                } else { throw ConfigurationError("total pool concentration is not a positive number"); }

                CodecMetric metric(value);
                metric.compile_barcode_tolerance(value, ontology);
            } else { throw ConfigurationError("codec element must be a dictionary"); }
        }

    }
};
bool Transcode::infer_ID(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
    buffer.clear();
    if(!decode_value_by_key< string >(key, buffer, container)) {
        if(infer_PU("PU", buffer, container, undetermined)) {
            encode_key_value(key, buffer, container, ontology);
            return true;
        } else { return false; }
    } else { return true; }
};
bool Transcode::infer_PU(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
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
void Transcode::compile_decoder_transformation(Value& value) {
    if(value.HasMember("transform")) {
        compile_transformation(value);

        Rule rule(decode_value_by_key< Rule >("transform", value));
        int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));

        /* validate all tokens refer to an existing input segment */
        for(auto& token : rule.token_array) {
            if(!(token.input_segment_index < input_segment_cardinality)) {
                throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
            }
            if(token.empty()) {
                throw ConfigurationError("token " + string(token) + " is empty");
            }
            if(!token.constant()) {
                throw ConfigurationError("token " + string(token) + " is not fixed width");
            }
        }

        /* annotate the decoder with cardinality information from the transfortmation */
        int32_t nucleotide_cardinality(0);
        vector< int32_t > barcode_length(rule.output_segment_cardinality, 0);
        for(auto& transform : rule.transform_array) {
            barcode_length[transform.output_segment_index] += transform.length();
            nucleotide_cardinality += transform.length();
        }
        encode_key_value("segment cardinality", rule.output_segment_cardinality, value, ontology);
        encode_key_value("nucleotide cardinality", nucleotide_cardinality, value, ontology);
        encode_key_value("barcode length", barcode_length, value, ontology);

        /* verify nucleotide cardinality and uniqueness
           and annotate each barcode element with the barcode segment cardinality */
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

        size_t segment_index(0);
        string barcode_segment;
        string barcode_sequence;
        set< string > unique_barcode_sequence;
        reference = value.FindMember("codec");
        if(reference != value.MemberEnd()) {
            if(reference->value.IsObject()) {
                Value& codec(reference->value);
                for(auto& record : codec.GetObject()) {
                    reference = record.value.FindMember("barcode");
                    if(reference !=  record.value.MemberEnd()) {
                        barcode_sequence.clear();
                        segment_index = 0;
                        for(const auto& element : reference->value.GetArray()) {
                            barcode_segment.assign(element.GetString(), element.GetStringLength());
                            if(static_cast< int32_t >(barcode_segment.size()) == barcode_length[segment_index]) {
                                barcode_sequence.append(barcode_segment);
                                ++segment_index;
                            } else {
                                string key(record.name.GetString(), record.name.GetStringLength());
                                string message;
                                message.append("expected ");
                                message.append(to_string(barcode_length[segment_index]));
                                message.append(" but found ");
                                message.append(to_string(barcode_segment.size()));
                                message.append(" nucleotides in segment ");
                                message.append(to_string(segment_index));
                                message.append(" of barcode ");
                                message.append(key);
                                throw ConfigurationError(message);
                            }
                        }
                        if(!unique_barcode_sequence.count(barcode_sequence)) {
                            unique_barcode_sequence.emplace(barcode_sequence);
                        } else {
                            throw ConfigurationError("duplicate barcode sequence " + barcode_sequence);
                        }
                    }
                    encode_key_value("segment cardinality", rule.output_segment_cardinality, record.value, ontology);
                }
            }
        }
        unique_barcode_sequence.clear();
    }
};
void Transcode::compile_output() {
    /* expand base output url path */
    standardize_url_value_by_key("base output url", ontology, ontology, IoDirection::OUT);
    URL base_output(decode_value_by_key< URL >("base output url", ontology));

    /* expand the report URL */
    standardize_url_value_by_key("report url", ontology, ontology, IoDirection::OUT);
    relocate_url_by_key("report url", ontology, ontology, base_output);

    /* load output template */
    compile_template();
    Rule rule(decode_value_by_key< Rule >("transform", ontology["template"]));

    /* encode the output_segment_cardinality deduced from the rule */
    const int32_t output_segment_cardinality(rule.output_segment_cardinality);
    encode_key_value("output segment cardinality", output_segment_cardinality, ontology, ontology);

    /* verify input segment references in the tokens reference existing input segments */
    const int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));
    for(auto& token : rule.token_array) {
        if(!(token.input_segment_index < input_segment_cardinality)) {
            throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
        }
    }

    Platform platform(decode_value_by_key< Platform >("platform", ontology));
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    uint8_t phred_offset(decode_value_by_key< uint8_t >("output phred offset", ontology));
    FormatType default_output_format(decode_value_by_key< FormatType >("default output format", ontology));
    FormatCompression default_output_compression(decode_value_by_key< FormatCompression >("default output compression", ontology));
    CompressionLevel default_output_compression_level(decode_value_by_key< CompressionLevel >("default output compression level", ontology));

    Value::MemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            Value& decoder_element(reference->value);

            /*  Collect channel references in a channel_reference_list
                This is a convience for iterating over channels in the
                codec element and the undetermined in once go */
            list< Value* > channel_reference_list;
            reference = decoder_element.FindMember("undetermined");
            if(reference != decoder_element.MemberEnd()) {
                channel_reference_list.push_back(&reference->value);
            }
            reference = decoder_element.FindMember("codec");
            if(reference != decoder_element.MemberEnd()) {
                for(auto& record : reference->value.GetObject()) {
                    channel_reference_list.push_back(&record.value);
                }
            }

            /* expand base output url path */
            standardize_url_value_by_key("base output url", decoder_element, ontology, IoDirection::OUT);
            URL base(decode_value_by_key< URL >("base output url", decoder_element));

            /* Collect query parameters from all reoccurring references of the same path */
            unordered_map< string, URL > canonical_url_by_path;
            for(auto& element : channel_reference_list) {
                standardize_url_array_by_key("output", *element, ontology, IoDirection::OUT);
                relocate_url_array_by_key("output", *element, ontology, base);
                list< URL > feed_url_array;
                if(decode_value_by_key< list< URL > >("output", feed_url_array, *element)) {
                    for(auto& url : feed_url_array) {
                        auto record = canonical_url_by_path.find(url.path());
                        if(record == canonical_url_by_path.end()) {
                            canonical_url_by_path[url.path()] = url;
                        } else {
                            record->second.override_query(url);
                        }
                    }
                }
            }

            /* Check for some obvious contradictions and set default type and compression */
            for(auto& url_record : canonical_url_by_path) {
                URL& url = url_record.second;
                if(url.is_stdin()) {
                    throw ConfigurationError("output stream can not be set to standard input");
                }
                if(url.is_stderr()) {
                    throw ConfigurationError("output stream can not be set to standard error");
                }
                if(url.type() == FormatType::UNKNOWN) {
                    url.set_type(default_output_format);
                }
                if(url.explicit_compression() == FormatCompression::UNKNOWN) {
                    url.set_compression(default_output_compression);
                }
                if(url.compression_level() == CompressionLevel::UNKNOWN) {
                    url.set_compression_level(default_output_compression_level);
                }
            }

            /* Apply the consensus query to all reoccurring references of the same path */
            for(auto& element : channel_reference_list) {
                list< URL > feed_url_array;
                if(decode_value_by_key< list< URL > >("output", feed_url_array, *element)) {
                    for(auto& url : feed_url_array) {
                        url = canonical_url_by_path[url.path()];
                    }
                }
                encode_key_value("output", feed_url_array, *element, ontology);
            }

            /* Assemble a map with the resolution of each feed in each channel */
            unordered_map< URL, unordered_map< int32_t, int > > feed_resolution;
            for(auto& element : channel_reference_list) {
                int32_t index(decode_value_by_key< int32_t >("index", *element));
                encode_key_value("TC", output_segment_cardinality, *element, ontology);
                pad_url_array_by_key("output", *element, output_segment_cardinality);
                list< URL > feed_url_array;
                if(decode_value_by_key< list< URL > >("output", feed_url_array, *element)) {
                    for(auto& url : feed_url_array) {
                        ++(feed_resolution[url][index]);
                    }
                }
            }
            if(feed_resolution.size() > 0) {
                int32_t index(0);
                unordered_map< URL, Value > feed_ontology_by_url;

                for(const auto& url_record : feed_resolution) {
                    const URL& url = url_record.first;

                    /* verify each feed url has the same resolution in all channels */
                    int resolution(0);
                    for(const auto& record : url_record.second) {
                        if(resolution == 0) {
                            resolution = record.second;
                        } else if(resolution != record.second) {
                            throw ConfigurationError("inconsistent resolution for " + string(url.path()));
                        }
                    }

                    /* make a proxy element for the feed */
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

                /* add the output feed by segment element to each channel */
                for(auto& element : channel_reference_list) {
                    list< URL > feed_url_array;
                    if(decode_value_by_key< list< URL > >("output", feed_url_array, *element)) {
                        Value feed_by_segment(kArrayType);
                        for(const auto& url : feed_url_array) {
                            const Value& proxy(feed_ontology_by_url[url]);
                            feed_by_segment.PushBack(Value(proxy, ontology.GetAllocator()).Move(), ontology.GetAllocator());
                        }
                        element->RemoveMember("feed");
                        element->AddMember("feed", Value(kObjectType).Move(), ontology.GetAllocator());
                        (*element)["feed"].AddMember("output feed by segment", feed_by_segment.Move(), ontology.GetAllocator());
                    }
                }

                Value feed_array(kArrayType);
                for(auto& record : feed_ontology_by_url) {
                    feed_array.PushBack(record.second.Move(), ontology.GetAllocator());
                }
                ontology["feed"].RemoveMember("output feed");
                ontology["feed"].AddMember("output feed", feed_array.Move(), ontology.GetAllocator());
            }
            cross_validate_io();
        }
    }
};
void Transcode::compile_template() {
    if(!ontology.HasMember("template")) {
        ontology.AddMember("template", Value(kObjectType).Move(), ontology.GetAllocator());
    }
    Value::MemberIterator reference = ontology.FindMember("template");
    Value& template_value(reference->value);

    if(!template_value.HasMember("transform")) {
        template_value.AddMember("transform", Value(kObjectType).Move(), ontology.GetAllocator());
    }

    reference = template_value.FindMember("transform");
    Value& transform_value(reference->value);

    const int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", ontology));

    /* if transform does not define a token array route all input segments to the output verbatim */
    if(!transform_value.HasMember("token")) {
        Value token_array(kArrayType);
        for(int32_t i(0); i < input_segment_cardinality; ++i) {
            string token(to_string(i) + "::");
            token_array.PushBack(Value(token.c_str(), token.size(), ontology.GetAllocator()), ontology.GetAllocator());
        }
        transform_value.AddMember("token", token_array.Move(), ontology.GetAllocator());
    }
    compile_transformation(template_value);
};
void Transcode::pad_url_array_by_key(const Value::Ch* key, Value& container, const int32_t& cardinality) {
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
void Transcode::cross_validate_io() {
    set< URL > input_url_set;
    Value::MemberIterator reference = ontology["feed"].FindMember("input feed");
    if(reference != ontology["feed"].MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            input_url_set.emplace(decode_value_by_key< URL >("url", element));
        }
    }
    set< URL > output_url_set;
    reference = ontology["feed"].FindMember("output feed");
    if(reference != ontology["feed"].MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            output_url_set.emplace(decode_value_by_key< URL >("url", element));
        }
    }
    URL report_url(decode_value_by_key< URL >("report url", ontology));
    if(!report_url.is_dev_null()) {
        if(input_url_set.count(report_url) > 0) {
            throw ConfigurationError("URL " + string(report_url) + " can not be used for both input and report");
        }
        if(output_url_set.count(report_url) > 0) {
            throw ConfigurationError("URL " + string(report_url) + " can not be used for both output and report");
        }
    }
    for(auto& url : output_url_set) {
        if(input_url_set.count(url) > 0) {
            throw ConfigurationError("URL " + string(url.path()) + " is used for both input and output");
        }
    }
};
void Transcode::compile_thread_model() {
    int32_t total_threads(decode_value_by_key< int32_t >("threads", ontology));

    int32_t decoding_threads;
    if(!decode_value_by_key("decoding threads", decoding_threads, ontology)) {
        decoding_threads = int32_t(round(double(total_threads) * (double(decoded_nucleotide_cardinality) / 1000.0)));
        decoding_threads = max(1, min(total_threads, max(1, decoding_threads)));
        encode_key_value("decoding threads", decoding_threads, ontology, ontology);
    }
};

/* validate */
void Transcode::validate() {
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
void Transcode::validate_decoder_group(const Value::Ch* key) {
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
void Transcode::validate_decoder(Value& value) {
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

/* execute */
void Transcode::load() {
    validate_url_accessibility();
    load_thread_pool();
    load_input();
    load_output();
    load_decoding();
    load_pivot();
};
void Transcode::validate_url_accessibility() {
    URL url;
    Value::MemberIterator reference = ontology["feed"].FindMember("input feed");
    if(reference != ontology["feed"].MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_readable()) {
                    throw IOError("can not open " + string(url.path()) + " for reading");
                }
            }
        }
    }

    reference = ontology["feed"].FindMember("output feed");
    if(reference != ontology["feed"].MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_writable()) {
                    throw IOError("can not open " + string(url.path()) + " for writing");
                }
            }
        }
    }
};
void Transcode::load_thread_pool() {
    if(thread_pool.pool == NULL) {
        int32_t htslib_threads(decode_value_by_key< int32_t >("htslib threads", ontology));
        thread_pool.pool = hts_tpool_init(htslib_threads);
        if(!thread_pool.pool) { throw InternalError("error creating thread pool"); }
    }
};
void Transcode::load_input() {
    if(input_feed_by_index.empty()) {
        /*  Decode feed_proxy_array, a local list of input feed proxy.
            The list has already been enumerated by the interface
            and contains only unique url references
        */
        list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("input feed", ontology["feed"]));

        /*  Initialized the hfile reference and verify input format */
        for(auto& proxy : feed_proxy_array) {
            proxy.open();
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
                throw InternalError("missing feed for URL " + string(url.path()) + " referenced in input proxy segment array");
            }
        }
    }
};
void Transcode::load_output() {
    /*  Decode feed_proxy_array, a local list of output feed proxy.
        The list has already been enumerated by the Pipeline
        and contains only unique url references
    */
    list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("output feed", ontology["feed"]));
    HeadPGAtom program(decode_value_by_key< HeadPGAtom >("program", ontology));

    /*  Register the read group elements on the feed proxy so it can be added to SAM header
        if a URL is present in the channel output that means the channel writes output to that file
        and the read group should be added to the header of that file.
    */
    map< URL, FeedProxy* > feed_proxy_by_url;
    for(auto& proxy : feed_proxy_array) {
        proxy.head.add_program(program);
        feed_proxy_by_url.emplace(make_pair(proxy.url, &proxy));
    };

    Value::ConstMemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                const Value& classifier(reference->value);

                reference = classifier.FindMember("undetermined");
                if(reference != classifier.MemberEnd()) {
                    HeadRGAtom rg(reference->value);
                    list< URL > output(decode_value_by_key< list< URL > >("output", reference->value));
                    for(auto& url : output) {
                        feed_proxy_by_url[url]->head.add_read_group(rg);
                    }
                }
                reference = classifier.FindMember("codec");
                if(reference != classifier.MemberEnd()) {
                    for(auto& record : reference->value.GetObject()) {
                        HeadRGAtom rg(record.value);
                        list< URL > output(decode_value_by_key< list< URL > >("output", record.value));
                        for(auto& url : output) {
                            feed_proxy_by_url[url]->head.add_read_group(rg);
                        }
                    }
                }
            } else { throw ConfigurationError("multiplex element must be a dictionary"); }
        }
    }

    /*  Initialized the hfile reference */
    for(auto& proxy : feed_proxy_array) {
        proxy.open();
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
                    throw ConfigurationError("unknown output format " + string(proxy.url));
                    break;
                };
            }
            feed->set_thread_pool(&thread_pool);
            output_feed_by_index.push_back(feed);
            output_feed_by_url.emplace(make_pair(proxy.url, feed));
    }
};
void Transcode::load_decoding() {
    load_multiplex_decoding();
    load_molecular_decoding();
    load_cellular_decoding();
};
void Transcode::load_multiplex_decoding() {
    Value::ConstMemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", reference->value));
        switch(algorithm) {
            case Algorithm::PAMLD: {
                sample_classifier = new PamlMultiplexDecoder(reference->value);
                break;
            };
            case Algorithm::MDD: {
                sample_classifier = new MdMultiplexDecoder(reference->value);
                break;
            };
            case Algorithm::TRANSPARENT: {
                sample_classifier = new RoutingClassifier< Channel >(reference->value);
                break;
            };
            default:
                throw ConfigurationError("unsupported multiplex decoder algorithm " + to_string(algorithm));
                break;
        }
    }
};
void Transcode::load_molecular_decoding() {
    Value::ConstMemberIterator reference = ontology.FindMember("molecular");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            molecular_classifier_array.reserve(1);
            load_molecular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            molecular_classifier_array.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_molecular_decoder(element);
            }
        }
    }
    molecular_classifier_array.shrink_to_fit();
};
void Transcode::load_molecular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch(algorithm) {
        case Algorithm::NAIVE: {
            molecular_classifier_array.emplace_back(new NaiveMolecularDecoder(value));
            break;
        };
        default:
            throw ConfigurationError("unsupported molecular decoder algorithm " + to_string(algorithm));
            break;
    }
};
void Transcode::load_cellular_decoding() {
    Value::ConstMemberIterator reference = ontology.FindMember("cellular");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            cellular_classifier_array.reserve(1);
            load_cellular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            cellular_classifier_array.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_cellular_decoder(element);
            }
        }
    }
    cellular_classifier_array.shrink_to_fit();
};
void Transcode::load_cellular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch(algorithm) {
        case Algorithm::PAMLD: {
            cellular_classifier_array.emplace_back(new PamlCellularDecoder(value));
            break;
        };
        case Algorithm::MDD: {
            cellular_classifier_array.emplace_back(new MdCellularDecoder(value));
            break;
        };
        default:
            throw ConfigurationError("unsupported cellular decoder algorithm " + to_string(algorithm));
            break;
    }
};
void Transcode::load_pivot() {
    int32_t decoding_threads(decode_value_by_key< int32_t >("decoding threads", ontology));
    for(int32_t index(0); index < decoding_threads; ++index) {
        pivot_array.emplace_back(*this, index);
    }
};
void Transcode::start() {
    for(auto feed : input_feed_by_index) {
        feed->initiate();
    }
    for(auto feed : output_feed_by_index) {
        feed->initiate();
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
void Transcode::stop() {
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
void Transcode::finalize() {
    Job::finalize();
    if(count > 0) {
        pf_fraction = double(pf_count) / double(count);
        Value element(kObjectType);
        encode_key_value("count", count, element, report);
        encode_key_value("pf count", pf_count, element, report);
        encode_key_value("pf fraction", pf_fraction, element, report);
        report.AddMember("incoming", element.Move(), report.GetAllocator());
    }

    /*  collect statistics from the accumulators on all pivot threads */
    for(auto& pivot : pivot_array) {
        *this += pivot;
    }

    if(sample_classifier != NULL) {
        sample_classifier->finalize();
        Value element(kObjectType);
        sample_classifier->encode(element, report);
        report.AddMember("multiplex", element.Move(), report.GetAllocator());
    }

    if(!molecular_classifier_array.empty()) {
        Value array(kArrayType);
        for(auto& classifier : molecular_classifier_array) {
            classifier->finalize();
            Value element(kObjectType);
            classifier->encode(element, report);
            array.PushBack(element.Move(), report.GetAllocator());
        }
        report.AddMember("molecular", array.Move(), report.GetAllocator());
    }

    if(!cellular_classifier_array.empty()) {
        Value array(kArrayType);
        for(auto& classifier : cellular_classifier_array) {
            classifier->finalize();
            Value element(kObjectType);
            classifier->encode(element, report);
            array.PushBack(element.Move(), report.GetAllocator());
        }
        report.AddMember("cellular", array.Move(), report.GetAllocator());
    }

    clean_json_value(report, report);
    sort_json_value(report, report);
};
void Transcode::apply_interactive_ontology(Document& document) const {
    Document adjusted;
    adjusted.CopyFrom(interactive, adjusted.GetAllocator());
    adjusted.RemoveMember("configuration url");
    adjusted.RemoveMember("static only");
    adjusted.RemoveMember("validate only");
    adjusted.RemoveMember("compile only");

    /* Format template token array into an output transform */
    Value::MemberIterator reference = adjusted.FindMember("template token");
    if(reference != adjusted.MemberEnd()) {
        Value template_value(kObjectType);
        Value transform_value(kObjectType);
        transform_value.AddMember("token", Value(reference->value, adjusted.GetAllocator()).Move(), adjusted.GetAllocator());
        template_value.AddMember("transform", transform_value.Move(), adjusted.GetAllocator());
        adjusted.AddMember("template", template_value.Move(), adjusted.GetAllocator());
    }
    adjusted.RemoveMember("template token");
    overlay_json_object(document, adjusted);
};

/* describe */
void Transcode::describe(ostream& o) const {
    print_global_instruction(o);
    print_input_instruction(o);
    print_transform_instruction(o);
    print_multiplex_instruction(o);
    print_molecular_instruction(o);
    print_cellular_instruction(o);
};
void Transcode::print_global_instruction(ostream& o) const {
    o << setprecision(float_precision());
    o << "Environment " << endl << endl;
    // o << "    Version                                     " << interface.application_version << endl;

    URL base_input_url(decode_value_by_key< URL >("base input url", ontology));
    o << "    Base input URL                              " << base_input_url << endl;

    URL base_output_url(decode_value_by_key< URL >("base input url", ontology));
    o << "    Base output URL                             " << base_output_url << endl;

    Platform platform(decode_value_by_key< Platform >("platform", ontology));
    o << "    Platform                                    " << platform << endl;

    bool enable_quality_control(decode_value_by_key< bool >("enable quality control", ontology));
    o << "    Quality tracking                            " << (enable_quality_control ? "enabled" : "disabled") << endl;

    bool filter_incoming_qc_fail(decode_value_by_key< bool >("filter incoming qc fail", ontology));
    o << "    Filter incoming QC failed reads             " << (filter_incoming_qc_fail ? "enabled" : "disabled") << endl;

    bool filter_outgoing_qc_fail(decode_value_by_key< bool >("filter outgoing qc fail", ontology));
    o << "    Filter outgoing QC failed reads             " << (filter_outgoing_qc_fail ? "enabled" : "disabled") << endl;

    uint8_t input_phred_offset(decode_value_by_key< uint8_t >("input phred offset", ontology));
    o << "    Input Phred offset                          " << to_string(input_phred_offset) << endl;

    uint8_t output_phred_offset(decode_value_by_key< uint8_t >("output phred offset", ontology));
    o << "    Output Phred offset                         " << to_string(output_phred_offset) << endl;

    int32_t leading_segment_index(decode_value_by_key< int32_t >("leading segment index", ontology));
    o << "    Leading segment index                       " << to_string(leading_segment_index) << endl;

    FormatType default_output_format(decode_value_by_key< FormatType >("default output format", ontology));
    o << "    Default output format                       " << default_output_format << endl;

    FormatCompression default_output_compression(decode_value_by_key< FormatCompression >("default output compression", ontology));
    o << "    Default output compression                  " << default_output_compression << endl;

    CompressionLevel default_output_compression_level(decode_value_by_key< CompressionLevel >("default output compression level", ontology));
    o << "    Default output compression level            " << default_output_compression_level << endl;

    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", ontology));
    o << "    Feed buffer capacity                        " << to_string(buffer_capacity) << endl;

    int32_t threads(decode_value_by_key< int32_t >("threads", ontology));
    o << "    Threads                                     " << to_string(threads) << endl;

    int32_t decoding_threads(decode_value_by_key< int32_t >("decoding threads", ontology));
    o << "    Decoding threads                            " << to_string(decoding_threads) << endl;

    int32_t htslib_threads(decode_value_by_key< int32_t >("htslib threads", ontology));
    o << "    HTSLib threads                              " << to_string(htslib_threads) << endl;
    o << endl;
};
void Transcode::print_feed_instruction(const Value::Ch* key, ostream& o) const {
    Value::ConstMemberIterator reference = ontology["feed"].FindMember(key);
    if(reference != ontology["feed"].MemberEnd()) {
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

                        switch (direction) {
                            case IoDirection::IN:
                                o << "    Input feed No." << index << endl;
                                o << "        Type : " << url.type() << endl;
                                if (url.compression() != FormatCompression::NONE) {
                                o << "        Compression : " << url.compression() << endl;
                                }
                                o << "        Resolution : " << resolution << endl;
                                o << "        Phred offset : " << to_string(phred_offset) << endl;
                                o << "        Platform : " << platform << endl;
                                o << "        Buffer capacity : " << capacity << endl;
                                o << "        URL : " << url << endl;
                                break;
                            case IoDirection::OUT:
                                o << "    Output feed No." << index << endl;
                                o << "        Type : " << url.type() << endl;
                                if (url.compression() != FormatCompression::NONE) {
                                o << "        Compression : " << url.compression() << "@" << url.compression_level() << endl;
                                }
                                o << "        Resolution : " << resolution << endl;
                                o << "        Phred offset : " << to_string(phred_offset) << endl;
                                o << "        Platform : " << platform << endl;
                                o << "        Buffer capacity : " << capacity << endl;
                                o << "        URL : " << url << endl;
                                break;
                            default:
                                break;
                        }
                        o << endl;
                    }
                }
            }
        }
    }
};
void Transcode::print_input_instruction(ostream& o) const {
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
void Transcode::print_transform_instruction(ostream& o) const {
    o << "Output transform" << endl << endl;

    int32_t output_segment_cardinality;
    if(decode_value_by_key< int32_t >("output segment cardinality", output_segment_cardinality, ontology)) {
        o << "    Output segment cardinality                  " << to_string(output_segment_cardinality) << endl;
    }

    Rule rule(decode_value_by_key< Rule >("transform", ontology["template"]));
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
void Transcode::print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const {
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
void Transcode::print_codec_instruction(const Value& value, const bool& plural, ostream& o) const {
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

            if(value.HasMember("transform")) {
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

                if(is_display_distance()) {
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
void Transcode::print_channel_instruction(const string& key, const Value& value, ostream& o) const {
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
            if(!barcode.empty()) { o << "        Barcode       : " << barcode << endl; }
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
void Transcode::print_multiplex_instruction(ostream& o) const {
    print_codec_group_instruction("multiplex", "Mutliplex decoding", o);
    print_feed_instruction("output feed", o);
};
void Transcode::print_molecular_instruction(ostream& o) const {
    print_codec_group_instruction("molecular", "Molecular decoding", o);
};
void Transcode::print_cellular_instruction(ostream& o) const {
    print_codec_group_instruction("cellular", "Cellular decoding", o);
};

TranscodePivot::TranscodePivot(Transcode& job, const int32_t& index) try :
    index(index),
    platform(decode_value_by_key< Platform >("platform", job.ontology)),
    leading_segment_index(decode_value_by_key< int32_t >("leading segment index", job.ontology)),
    input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", job.ontology)),
    output_segment_cardinality(decode_value_by_key< int32_t >("output segment cardinality", job.ontology)),
    input(input_segment_cardinality, platform, leading_segment_index),
    output(output_segment_cardinality, platform, leading_segment_index),
    sample_classifier(NULL),
    job(job),
    filter_incoming_qc_fail(decode_value_by_key< bool >("filter incoming qc fail", job.ontology)),
    enable_quality_control(decode_value_by_key< bool >("enable quality control", job.ontology)),
    template_rule(decode_value_by_key< Rule >("transform", job.ontology["template"])) {

    load_decoding();

    } catch(Error& error) {
        error.push("TranscodePivot");
        throw;
};
void TranscodePivot::load_decoding() {
    load_multiplex_decoding();
    load_molecular_decoding();
    load_cellular_decoding();
    input.clear();
    output.clear();
};
void TranscodePivot::load_multiplex_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("multiplex");
    if(reference != job.ontology.MemberEnd()) {
        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", reference->value));
        switch (algorithm) {
            case Algorithm::PAMLD: {
                PamlMultiplexDecoder* pamld_decoder(new PamlMultiplexDecoder(reference->value));
                pamld_decoder->unclassified.populate(job.output_feed_by_url);
                for(auto& channel : pamld_decoder->element_by_index) {
                    channel.populate(job.output_feed_by_url);
                }
                sample_classifier = pamld_decoder;
                break;
            };
            case Algorithm::MDD: {
                MdMultiplexDecoder* mdd_decoder(new MdMultiplexDecoder(reference->value));
                mdd_decoder->unclassified.populate(job.output_feed_by_url);
                for(auto& channel : mdd_decoder->element_by_index) {
                    channel.populate(job.output_feed_by_url);
                }
                sample_classifier = mdd_decoder;
                break;
            };
            case Algorithm::TRANSPARENT: {
                RoutingClassifier< Channel >* transparent_decoder(new RoutingClassifier< Channel >(reference->value));
                transparent_decoder->unclassified.populate(job.output_feed_by_url);
                sample_classifier = transparent_decoder;
                break;
            };
            default:
                throw ConfigurationError("unknown multiplex decoder algorithm");
                break;
        }
    }
};
void TranscodePivot::load_molecular_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("molecular");
    if(reference != job.ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            molecular_classifier_array.reserve(1);
            load_molecular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            molecular_classifier_array.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_molecular_decoder(element);
            }
        }
    }
    molecular_classifier_array.shrink_to_fit();
};
void TranscodePivot::load_molecular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch (algorithm) {
        case Algorithm::NAIVE: {
            NaiveMolecularDecoder* naive_decoder(new NaiveMolecularDecoder(value));
            molecular_classifier_array.emplace_back(naive_decoder);
            break;
        };
        default:
            break;
    }
};
void TranscodePivot::load_cellular_decoding() {
    Value::ConstMemberIterator reference = job.ontology.FindMember("cellular");
    if(reference != job.ontology.MemberEnd()) {
        if(reference->value.IsObject()) {
            cellular_classifier_array.reserve(1);
            load_cellular_decoder(reference->value);

        } else if(reference->value.IsArray()) {
            cellular_classifier_array.reserve(reference->value.Size());
            for(const auto& element : reference->value.GetArray()) {
                load_cellular_decoder(element);
            }
        }
    }
    cellular_classifier_array.shrink_to_fit();
};
void TranscodePivot::load_cellular_decoder(const Value& value) {
    Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
    switch (algorithm) {
        case Algorithm::PAMLD: {
            PamlCellularDecoder* paml_decoder(new PamlCellularDecoder(value));
            cellular_classifier_array.emplace_back(paml_decoder);
            break;
        };
        case Algorithm::MDD: {
            MdCellularDecoder* md_decoder(new MdCellularDecoder(value));
            cellular_classifier_array.emplace_back(md_decoder);
            break;
        };
        default:
            break;
    }
};
