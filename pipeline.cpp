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

Pivot::Pivot(Pipeline& pipeline, const Value& ontology, const int32_t& index) :
    index(index),
    input(pipeline.input_segment_cardinality, pipeline.platform, pipeline.leading_segment_index),
    output(pipeline.output_segment_cardinality, pipeline.platform, pipeline.leading_segment_index),
    multiplex(NULL),
    input_accumulator(ontology),
    output_accumulator(find_value_by_key("multiplex", ontology)),
    pipeline(pipeline),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", ontology)),
    template_rule(decode_value_by_key< Rule >("template", ontology)) {

    load_multiplex_decoder(ontology);
    load_molecular_decoder(ontology);
    load_splitseq_decoder(ontology);
    clear();
};
void Pivot::load_multiplex_decoder(const Value& ontology) {
    Value::ConstMemberIterator reference = ontology.FindMember("multiplex");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", reference->value));
                BarcodeDecoder< Channel >* decoder;
                switch (algorithm) {
                    case Algorithm::PAMLD: {
                        decoder = new PAMLMultiplexDecoder(reference->value);
                        break;
                    };
                    case Algorithm::MDD: {
                        decoder = new MDMultiplexDecoder(reference->value);
                        break;
                    };
                    default:
                        throw ConfigurationError("unknown multiplex decoder algorithm");
                        break;
                }
                pipeline.populate_multiplex_decoder(*decoder);
                multiplex = decoder;
            } else { throw ConfigurationError("multiplex element must be a dictionary"); }
        }
    }
};
void Pivot::load_molecular_decoder(const Value& ontology) {
    Value::ConstMemberIterator reference = ontology.FindMember("molecular");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                molecular.reserve(reference->value.Size());
                for(const auto& element : reference->value.GetArray()) {
                    if(element.IsObject()) {
                        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", element));
                        switch (algorithm) {
                            case Algorithm::SIMPLE: {
                                molecular.emplace_back(new SimpleMolecularDecoder(element));
                                break;
                            };
                            default:
                                throw ConfigurationError("unknown molecular decoder algorithm");
                                break;
                        }
                    } else { throw ConfigurationError("molecular decoder array element must be a dictionary"); }


                }
            } else { throw ConfigurationError("molecular decoder element must be an array"); }
        }
    }
};
void Pivot::load_splitseq_decoder(const Value& ontology) {
    Value::ConstMemberIterator reference = ontology.FindMember("splitseq");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                splitseq.reserve(reference->value.Size());
                for(const auto& element : reference->value.GetArray()) {
                    if(element.IsObject()) {
                        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", element));
                        switch (algorithm) {
                            case Algorithm::PAMLD: {
                                splitseq.emplace_back(new PAMLSplitSEQDecoder(element));
                                break;
                            };
                            case Algorithm::MDD: {
                                splitseq.emplace_back(new MDSplitSEQDecoder(element));
                                break;
                            };
                            default:
                                throw ConfigurationError("unknown SplitSEQ decoder algorithm");
                                break;
                        }
                    } else { throw ConfigurationError("SplitSEQ decoder array element must be a dictionary"); }
                }
            } else { throw ConfigurationError("SplitSEQ decoder element must be an array"); }
        }
    }
};
void Pivot::run() {
    while(pipeline.pull(input)) {
        input.validate();
        transform();
        push();
        increment();
        clear();
    }
};

/*  Pipeline */

Pipeline::Pipeline(const Document& instruction) :
    instruction(instruction),
    platform(decode_value_by_key< Platform >("platform", instruction)),
    leading_segment_index(decode_value_by_key< int32_t >("leading segment index", instruction)),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", instruction)),
    include_filtered(decode_value_by_key< bool >("include filtered", instruction)),
    threads(decode_value_by_key< int32_t >("threads", instruction)),
    input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", instruction)),
    output_segment_cardinality(decode_value_by_key< int32_t >("output segment cardinality", instruction)),
    end_of_input(false),
    thread_pool({NULL, 0}),
    input_accumulator(instruction),
    output_accumulator(find_value_by_key("multiplex", instruction)) {

    thread_pool.pool = hts_tpool_init(threads);
    if(!thread_pool.pool) { throw InternalError("error creating thread pool"); }
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
};
bool Pipeline::pull(Read& read) {
    vector< unique_lock< mutex > > feed_locks;
    feed_locks.reserve(input_feed_by_index.size());

    // acquire a pull lock for all feeds in a fixed order
    for(const auto feed : input_feed_by_index) {
        feed_locks.push_back(feed->acquire_pull_lock());
    }

    // pull into pivot input segments from input feeds
    for(size_t i = 0; i < read.segment_cardinality(); ++i) {
        if(!input_feed_by_segment[i]->pull(read[i])) {
            end_of_input = true;
        }
    }

    // release the locks on the feeds in reverse order
    for(auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
        feed_lock->unlock();
    }
    return !end_of_input;
};
void Pipeline::execute() {
    load();
    start();
    stop();
    finalize();
};
void Pipeline::load() {
    validate_url_accessibility();
    load_input();
    load_output();
    load_pivot();
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
    for(auto& pivot : pivot_array) {
        pivot.start();
    }
    for(auto& pivot : pivot_array) {
        pivot.join();
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
    for(auto& pivot : pivot_array) {
        input_accumulator += pivot.input_accumulator;
        output_accumulator += pivot.output_accumulator;
    }
    input_accumulator.finalize();
    output_accumulator.finalize();
    encode_report(cout);
};
void Pipeline::validate_url_accessibility() {
    URL url;
    Value::ConstMemberIterator reference = instruction.FindMember("input feed");
    if(reference != instruction.MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_readable()) {
                    throw IOError("could not open " + string(url) + " for reading");
                }
            }
        }
    }

    reference = instruction.FindMember("output feed");
    if(reference != instruction.MemberEnd()) {
        for(auto& element : reference->value.GetArray()) {
            if(decode_value_by_key< URL >("url", url, element)) {
                if(!url.is_writable()) {
                    throw IOError("could not open " + string(url) + " for writing");
                }
            }
        }
    }
};
void Pipeline::load_input() {
    /*  Decode feed_proxy_array, a local list of input feed proxy.
        The list has already been enumerated by the interface
        and contains only unique url references
    */
    list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("input feed", instruction));

    /*  Initialized the hfile reference and verify input format */
    for(auto& proxy : feed_proxy_array) {
        proxy.probe();
    };

    /*  Load feed_by_url, a local map of input feeds by url, from the proxy.
        Populate input_feed_by_index used to enumerate threaded access to the input feeds */
    unordered_map< URL, Feed* > feed_by_url(feed_proxy_array.size());
    for(auto& proxy : feed_proxy_array) {
        Feed* feed(NULL);
        FormatKind kind(format_kind_from_type(proxy.url.type()));
        switch(kind) {
            case FormatKind::FASTQ: {
                feed = new FastqFeed(proxy);
                break;
            };
            case FormatKind::HTS: {
                feed = new HtsFeed(proxy);
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

    list< URL > url_by_segment(decode_value_by_key< list < URL > >("input", instruction));

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
};
void Pipeline::load_output() {
    /*  Decode feed_proxy_array, a local list of output feed proxy.
        The list has already been enumerated by the environment
        and contains only unique url references
    */
    list< FeedProxy > feed_proxy_array(decode_value_by_key< list< FeedProxy > >("output feed", instruction));

    /*  Register the read group elements on the feed proxy so it can be added to SAM header
        if a URL is present in the channel output that means the channel writes output to that file
        and the read group should be added to the header of that file.
    */
    map< URL, FeedProxy* > feed_proxy_by_url;
    for(auto& proxy : feed_proxy_array) {
        feed_proxy_by_url.emplace(make_pair(proxy.url, &proxy));
    };

    Value::ConstMemberIterator reference = instruction.FindMember("multiplex");
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                const Value& multiplex(reference->value);
                reference = multiplex.FindMember("codec");
                if(reference != instruction.MemberEnd()) {
                    if(!reference->value.IsNull()) {
                        if(reference->value.IsObject()) {
                            for(auto& record : reference->value.GetObject()) {
                                HeadRGAtom rg(record.value);
                                list< URL > output(decode_value_by_key< list< URL > >("output", record.value));
                                for(auto& url : output) {
                                    feed_proxy_by_url[url]->register_rg(rg);
                                }
                            }
                        } else { throw ConfigurationError("codec element must be a dictionary"); }
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
        if(!proxy.url.is_null()) {
            Feed* feed(NULL);
            FormatKind kind(format_kind_from_type(proxy.url.type()));
            switch(kind) {
                case FormatKind::FASTQ: {
                    feed = new FastqFeed(proxy);
                    break;
                };
                case FormatKind::HTS: {
                    feed = new HtsFeed(proxy);
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
    }
};
void Pipeline::load_pivot() {
    for(int32_t index(0); index < threads; ++index) {
        pivot_array.emplace_back(*this, instruction, index);
    }
};
void Pipeline::populate_multiplex_decoder(DiscreteDecoder< Channel >& decoder) {
    populate_multiplex_channel(decoder.undetermined);
    for(auto& channel : decoder.element_by_index) {
        populate_multiplex_channel(channel);
    }
};
void Pipeline::populate_multiplex_channel(Channel& channel) {
    map< int32_t, Feed* > feed_by_index;
    channel.output_feed_by_segment.reserve(channel.output_feed_url_by_segment.size());
    for(const auto& url : channel.output_feed_url_by_segment) {
        if(!url.is_null()) {
            Feed* feed(output_feed_by_url[url]);
            channel.output_feed_by_segment.emplace_back(feed);
            if(feed_by_index.count(feed->index) == 0) {
                feed_by_index.emplace(make_pair(feed->index, feed));
            }
        }
    }
    channel.output_feed_by_segment.shrink_to_fit();

    channel.output_feed_by_order.reserve(feed_by_index.size());
    for(auto& record : feed_by_index) {
        channel.output_feed_by_order.push_back(record.second);
    }
    channel.output_feed_by_order.shrink_to_fit();
};
void Pipeline::encode_report(ostream& o) const {
    Document document(kObjectType);

    encode_key_value("demultiplex output report", output_accumulator, document, document);
    encode_key_value("demultiplex input report", input_accumulator, document, document);

    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    document.Accept(writer);
    o << buffer.GetString();
    o << endl;
};
