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

#include "interface.h"

void to_string(const ProgramAction& value, string& result) {
    switch(value) {
        case ProgramAction::DEMULTIPLEX:    result.assign("demux");      break;
        case ProgramAction::QUALITY:        result.assign("quality");    break;
        default:                            result.assign("unknown");    break;
    }
};
string to_string(const ProgramAction& value) {
    string string_value;
    to_string(value, string_value);
    return string_value;
};
bool from_string(const char* value, ProgramAction& result) {
         if(value == NULL)              result = ProgramAction::UNKNOWN;
    else if(!strcmp(value, "demux"))    result = ProgramAction::DEMULTIPLEX;
    else if(!strcmp(value, "quality"))  result = ProgramAction::QUALITY;
    else                                result = ProgramAction::UNKNOWN;
    return (result == ProgramAction::UNKNOWN ? false : true);
};
bool from_string(const string& value, ProgramAction& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const ProgramAction& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const ProgramAction& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< ProgramAction >(const Value::Ch* key, ProgramAction& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

static int32_t inheritence_depth(const string& key, const unordered_map< string, Value* >& node_by_key, Document& document) {
    int32_t depth(0);
    auto record = node_by_key.find(key);
    if(record != node_by_key.end()) {
        Value* value = record->second;
        if(!decode_value_by_key("depth", depth, *value)) {
            string base_key;
            if(decode_value_by_key("base", base_key, *value)) {
                if(base_key != key) {
                    depth = inheritence_depth(base_key, node_by_key, document) + 1;
                    encode_key_value("depth", depth, *value, document);
                } else { throw ConfigurationError("object can not inherit from itself " + key); }
            } else {
                encode_key_value("depth", depth, *value, document);
            }
        }
    } else { throw ConfigurationError("referencing an undefined base decoder " + key); }
    return depth;
};

void Demultiplex::compile_instruction() {
    compile_PG_element();
    compile_input_instruction();
    apply_decoder_inheritence();
    compile_decoder_group_instruction("multiplex");
    compile_decoder_group_instruction("molecular");
    compile_decoder_group_instruction("splitseq");
    compile_output_instruction();
};
void Demultiplex::compile_transformation(Value& value) {
    /* add the default observation if one was not specificed.
       default observation will treat every token as a segment */
    if(value.IsObject()) {
        Value::MemberIterator reference = value.FindMember("template");
        if(reference != value.MemberEnd()) {
            if(reference->value.IsObject()) {
                Value& template_element(reference->value);
                reference = template_element.FindMember("token");
                if(reference != template_element.MemberEnd()) {
                    if(reference->value.IsArray()) {
                        int32_t token_cardinality(reference->value.Size());

                        reference = template_element.FindMember("observation");
                        if(reference == template_element.MemberEnd() || reference->value.IsNull() || (reference->value.IsArray() && reference->value.Empty())) {
                            Value observation(kArrayType);
                            for(int32_t i(0); i < token_cardinality; ++i) {
                                string element(to_string(i));
                                observation.PushBack(Value(element.c_str(), element.size(),instruction.GetAllocator()).Move(), instruction.GetAllocator());
                            }
                            template_element.RemoveMember("observation");
                            template_element.AddMember(Value("observation", instruction.GetAllocator()).Move(), observation.Move(), instruction.GetAllocator());
                        }
                    } else { throw ConfigurationError("template token element is not an array"); }
                } else { throw ConfigurationError("template element is missing a token array"); }
            }
        }
    }
};
void Demultiplex::apply_decoder_inheritence() {
    Value::MemberIterator reference = instruction.FindMember("decoder");
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {

                /*  map each decoder by key */
                unordered_map< string, Value* > codec_by_key;
                for(auto& record : reference->value.GetObject()) {
                    if(!record.value.IsNull()) {
                        record.value.RemoveMember("depth");
                        codec_by_key.emplace(make_pair(string(record.name.GetString(), record.name.GetStringLength()), &record.value));
                    }
                }

                /* compute the inheritence depth of each decoder and keep track of the max depth */
                int32_t max_depth(0);
                for(auto& record : codec_by_key) {
                    try {
                        max_depth = max(max_depth, inheritence_depth(record.first, codec_by_key, instruction));
                    } catch(ConfigurationError& error) {
                        throw CommandLineError("decoder " + record.first + " is " + error.message);
                    }
                }

                /* apply decoder inheritence down the tree */
                int32_t depth(0);
                for(int32_t i(1); i <= max_depth; ++i) {
                    for(auto& record : codec_by_key) {
                        Value* value = record.second;
                        if(decode_value_by_key("depth", depth, *value) && depth == i) {
                            string base;
                            if(decode_value_by_key< string >("base", base, *value)) {
                                merge_json_value(*codec_by_key[base], *value, instruction);
                            }
                        }
                    }
                }
            } else { throw ConfigurationError("decoder element must be a dictionary"); }
        }
    }
};
void Demultiplex::compile_input_instruction() {
    expand_url_value_by_key("base input url", instruction, instruction);
    expand_url_array_by_key("input", instruction, instruction, IoDirection::IN);
    URL base(decode_value_by_key< URL >("base input url", instruction));

    relocate_url_array_by_key("input", instruction, instruction, base);
    list< URL > feed_url_array(decode_value_by_key< list< URL > >("input", instruction));

    Platform platform(decode_value_by_key< Platform >("platform", instruction));
    int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", instruction));
    uint8_t input_phred_offset(decode_value_by_key< uint8_t >("input phred offset", instruction));

    /* encode the input segment cardinality */
    int32_t input_segment_cardinality(static_cast< int32_t>(feed_url_array.size()));
    encode_key_value("input segment cardinality", input_segment_cardinality, instruction, instruction);

    /*  validate leading_segment_index */
    int32_t leading_segment_index(decode_value_by_key< int32_t >("leading segment index", leading_segment_index, instruction));
    if(leading_segment_index >= input_segment_cardinality) {
        throw ConfigurationError("leading segment index " + to_string(leading_segment_index) + " references non existing input segment");
    }

    map< URL, int > feed_resolution;
    for(const auto& url : feed_url_array) {
        ++(feed_resolution[url]);
    }

    int32_t feed_index(0);
    unordered_map< URL, Value > feed_ontology_by_url;
    for(const auto& record : feed_resolution) {
        const URL& url = record.first;
        int resolution(record.second);
        Value proxy(kObjectType);
        encode_key_value("index", feed_index, proxy, instruction);
        encode_key_value("url", record.first, proxy, instruction);
        encode_key_value("direction", IoDirection::IN, proxy, instruction);
        encode_key_value("platform", platform, proxy, instruction);
        encode_key_value("capacity", buffer_capacity * resolution, proxy, instruction);
        encode_key_value("resolution", resolution, proxy, instruction);
        encode_key_value("phred offset", input_phred_offset, proxy, instruction);
        feed_ontology_by_url.emplace(make_pair(url, move(proxy)));
        ++feed_index;
    }

    Value feed_by_segment(kArrayType);
    for(const auto& url : feed_url_array) {
        const Value& proxy(feed_ontology_by_url[url]);
        feed_by_segment.PushBack(Value(proxy, instruction.GetAllocator()).Move(), instruction.GetAllocator());
    }
    instruction.RemoveMember("input feed by segment");
    instruction.AddMember("input feed by segment", feed_by_segment.Move(), instruction.GetAllocator());

    Value feed_array(kArrayType);
    for(auto& record : feed_ontology_by_url) {
        feed_array.PushBack(record.second.Move(), instruction.GetAllocator());
    }
    instruction.RemoveMember("input feed");
    instruction.AddMember("input feed", feed_array.Move(), instruction.GetAllocator());
};
void Demultiplex::compile_decoder_group_instruction(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            Value& node(reference->value);

            /* get the base decoder and codec node from the configuration */
            Value default_configuration_decoder(kObjectType);
            Value default_configuration_codec(kObjectType);
            Value::ConstMemberIterator const_reference = ontology.FindMember("projection");
            if(const_reference != ontology.MemberEnd()) {
                const Value& base = const_reference->value;
                if(!base.IsNull()) {
                    const_reference = base.FindMember("decoder");
                    if(const_reference != base.MemberEnd() && !const_reference->value.IsNull()) {
                        default_configuration_decoder.CopyFrom(const_reference->value, instruction.GetAllocator());
                    }
                    const_reference = base.FindMember("codec");
                    if(const_reference != base.MemberEnd() && !const_reference->value.IsNull()) {
                        default_configuration_codec.CopyFrom(const_reference->value, instruction.GetAllocator());
                    }
                }
            }

            /* project decoder attributes from the document root */
            Value default_instruction_decoder(kObjectType);
            project_json_value(default_configuration_decoder, instruction, default_instruction_decoder, instruction);
            default_configuration_decoder.SetNull();

            /* project codec attributes from the instruction root */
            Value default_instruction_codec(kObjectType);
            project_json_value(default_configuration_codec, instruction, default_instruction_codec, instruction);
            default_configuration_codec.SetNull();

            /* create a map of globally available codecs */
            unordered_map< string, Value* > codec_by_key;
            reference = instruction.FindMember("decoder");
            if(reference != instruction.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    if(reference->value.IsObject()) {
                        for(auto& record : reference->value.GetObject()) {
                            if(!record.value.IsNull()) {
                                string key(record.name.GetString(), record.name.GetStringLength());
                                codec_by_key.emplace(make_pair(key, &record.value));
                            }
                        }
                    }
                }
            }

            if(node.IsObject()) {
                string base;
                if(decode_value_by_key< string >("base", base, node)) {
                    auto record = codec_by_key.find(base);
                    if(record != codec_by_key.end()) {
                        merge_json_value(*record->second, node, instruction);
                    } else { throw ConfigurationError(string(key) + " is referencing an undefined base decoder " + base); }
                }
                encode_key_value("index", 0, node, instruction);
                compile_decoder_codec(node, default_instruction_decoder, default_instruction_codec);
                compile_decoder_transformation(node);

            } else if(node.IsArray()) {
                int32_t index(0);
                for(auto& element : node.GetArray()) {
                    if(element.IsObject()) {
                        string base;
                        if(decode_value_by_key< string >("base", base, element)) {
                            auto record = codec_by_key.find(base);
                            if(record != codec_by_key.end()) {
                                merge_json_value(*record->second, element, instruction);
                            } else { throw ConfigurationError("element " + to_string(index) + " of " + string(key) + " is referencing an undefined base decoder " + base); }
                        }
                        encode_key_value("index", index, element, instruction);
                        compile_decoder_codec(element, default_instruction_decoder, default_instruction_codec);
                        compile_decoder_transformation(element);
                        ++index;
                    } else { throw ConfigurationError("decoder element must be a dictionary"); }
                }
            }
            // compile_decoder_group(key);

            default_instruction_decoder.SetNull();
            default_instruction_codec.SetNull();
        }
    }
};
void Demultiplex::compile_decoder_codec(Value& value, const Value& default_instruction_decoder, const Value& default_instruction_codec) {
    if(value.IsObject()) {
        string buffer;

        /* merge the default instruction decoder */
        merge_json_value(default_instruction_decoder, value, instruction);

        /* compute default barcode induced by the codec */
        Value default_codec(kObjectType);
        project_json_value(default_instruction_codec, value, default_codec, instruction);
        clean_json_value(value, instruction);


        int32_t barcode_index(0);
        double total_concentration(0);
        set< string > unique_barcode_id;
        double noise(decode_value_by_key< double >("noise", value));

        /* apply barcode default on undetermined barcode or create it from the default if one was not explicitly specified */
        Value::MemberIterator reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()){
            merge_json_value(default_codec, reference->value, instruction);
        } else {
            value.AddMember (
                Value("undetermined", instruction.GetAllocator()).Move(),
                Value(default_codec, instruction.GetAllocator()).Move(),
                instruction.GetAllocator()
            );
        }
        clean_json_value(reference->value, instruction);

        reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()){
            encode_key_value("index", barcode_index, reference->value, instruction);
            if(infer_ID("ID", buffer, reference->value, true)) {
                unique_barcode_id.emplace(buffer);
            }
            encode_key_value("concentration", noise, reference->value, instruction);
            ++barcode_index;
        }

        reference = value.FindMember("codec");
        if(reference != value.MemberEnd()){
            if(reference->value.IsObject()) {
                Value& codec(reference->value);

                for(auto& record : codec.GetObject()) {
                    merge_json_value(default_codec, record.value, instruction);
                    clean_json_value(record.value, instruction);

                    encode_key_value("index", barcode_index, record.value, instruction);

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
                        encode_key_value("concentration", concentration * factor, record.value, instruction);
                    }
                } else { throw ConfigurationError("total pool concentration is not a positive number"); }
            } else { throw ConfigurationError("codec element must be a dictionary"); }
        }
        default_codec.SetNull();
    }
};
void Demultiplex::compile_decoder_transformation(Value& value) {
    compile_transformation(value);

    /* decode the transformation rule */
    Rule rule(decode_value_by_key< Rule >("template", value));
    int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", instruction));

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
    encode_key_value("segment cardinality", rule.output_segment_cardinality, value, instruction);
    encode_key_value("nucleotide cardinality", nucleotide_cardinality, value, instruction);
    encode_key_value("barcode length", barcode_length, value, instruction);

    /* annotate each barcode element with the barcode segment cardinality */
    Value::MemberIterator reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            Value& undetermined(reference->value);

            /* explicitly define a null barcode segment for the right dimension in the undetermined */
            Value barcode_segment(kArrayType);
            for(size_t i = 0; i < barcode_length.size(); ++i) {
                string sequence(barcode_length[i], '=');
                barcode_segment.PushBack(Value(sequence.c_str(), sequence.size(), instruction.GetAllocator()).Move(), instruction.GetAllocator());
            }
            undetermined.RemoveMember("barcode segment");
            undetermined.AddMember(Value("barcode segment", instruction.GetAllocator()).Move(), barcode_segment.Move(), instruction.GetAllocator());
            encode_key_value("segment cardinality", rule.output_segment_cardinality, undetermined, instruction);
        }
    }

    reference = value.FindMember("codec");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                for(auto& record : reference->value.GetObject()) {
                    encode_key_value("segment cardinality", rule.output_segment_cardinality, record.value, instruction);
                }
            }
        }
    }
};
bool Demultiplex::infer_PU(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
    string suffix;
    if(!decode_value_by_key< string >(key, suffix, container)) {
        if(!undetermined) {
            list< string > barcode;
            if(decode_value_by_key< list< string > >("barcode segment", barcode, container)) {
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
            encode_key_value(key, buffer, container, instruction);
            return true;
        } else { return false; }
    } else { return true; }
};
bool Demultiplex::infer_ID(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined) {
    if(!decode_value_by_key< string >(key, buffer, container)) {
        if(infer_PU("PU", buffer, container, undetermined)) {
            encode_key_value(key, buffer, container, instruction);
            return true;
        } else { return false; }
    } else { return true; }
};
void Demultiplex::pad_url_array_by_key(const Value::Ch* key, Value& container, const int32_t& cardinality) {
    list< URL > array;
    if(decode_value_by_key< list< URL > >(key, array, container)) {
        if(!array.empty()) {
            if(static_cast< int32_t >(array.size()) != cardinality) {
                if(array.size() == 1) {
                    while(static_cast< int32_t >(array.size()) < cardinality) {
                        array.push_back(array.front());
                    }
                    encode_key_value(key, array, container, instruction);
                } else { throw ConfigurationError("incorrect number of output URLs in channel"); }
            }
        }
    }
};
void Demultiplex::compile_output_instruction() {
    /* load output template */
    const int32_t input_segment_cardinality(decode_value_by_key< int32_t >("input segment cardinality", instruction));
    compile_transformation(instruction);

    Rule rule(decode_value_by_key< Rule >("template", instruction));
    for(auto& token : rule.token_array) {
        if(!(token.input_segment_index < input_segment_cardinality)) {
            throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
        }
    }
    const int32_t output_segment_cardinality(rule.output_segment_cardinality);
    encode_key_value("output segment cardinality", output_segment_cardinality, instruction, instruction);

    Value::MemberIterator reference = instruction.FindMember("multiplex");
    if(reference != instruction.MemberEnd()) {
        if(reference->value.IsObject()) {
            Value& value(reference->value);
            expand_url_value_by_key("base output url", value, instruction);
            URL base(decode_value_by_key< URL >("base output url", value));

            Platform platform(decode_value_by_key< Platform >("platform", instruction));
            int32_t buffer_capacity(decode_value_by_key< int32_t >("buffer capacity", instruction));
            uint8_t phred_offset(decode_value_by_key< uint8_t >("output phred offset", instruction));

            unordered_map< URL, unordered_map< int32_t, int > > feed_resolution;
            Value::MemberIterator reference = value.FindMember("undetermined");
            if(reference != value.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    int32_t index(decode_value_by_key< int32_t >("index", reference->value));
                    encode_key_value("TC", output_segment_cardinality, reference->value, instruction);
                    pad_url_array_by_key("output", reference->value, output_segment_cardinality);
                    expand_url_array_by_key("output", reference->value, instruction, IoDirection::OUT);
                    relocate_url_array_by_key("output", reference->value, instruction, base);
                    list< URL > feed_url_array(decode_value_by_key< list< URL > >("output", reference->value));

                    for(auto& url : feed_url_array) {
                        ++(feed_resolution[url][index]);
                    }
                }
            }
            reference = value.FindMember("codec");
            if(reference != value.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    if(reference->value.IsObject()) {
                        for(auto& record : reference->value.GetObject()) {
                            int32_t index(decode_value_by_key< int32_t >("index", record.value));
                            encode_key_value("TC", output_segment_cardinality, record.value, instruction);
                            pad_url_array_by_key("output", record.value, output_segment_cardinality);
                            expand_url_array_by_key("output", record.value, instruction, IoDirection::OUT);
                            relocate_url_array_by_key("output", record.value, instruction, base);
                            list< URL > feed_url_array(decode_value_by_key< list< URL > >("output", record.value));

                            for(auto& url : feed_url_array) {
                                ++(feed_resolution[url][index]);
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
                    encode_key_value("index", index, proxy, instruction);
                    encode_key_value("url", url, proxy, instruction);
                    encode_key_value("direction", IoDirection::OUT, proxy, instruction);
                    encode_key_value("platform", platform, proxy, instruction);
                    encode_key_value("capacity", buffer_capacity * resolution, proxy, instruction);
                    encode_key_value("resolution", resolution, proxy, instruction);
                    encode_key_value("phred offset", phred_offset, proxy, instruction);
                    feed_ontology_by_url.emplace(make_pair(url, move(proxy)));
                    ++index;
                }

                reference = value.FindMember("undetermined");
                if(reference != value.MemberEnd()) {
                    if(!reference->value.IsNull()) {
                        if(reference->value.IsObject()) {
                            list< URL > feed_url_array(decode_value_by_key< list< URL > >("output", reference->value));
                            Value feed_by_segment(kArrayType);
                            for(const auto& url : feed_url_array) {
                                const Value& proxy(feed_ontology_by_url[url]);
                                feed_by_segment.PushBack(Value(proxy, instruction.GetAllocator()).Move(), instruction.GetAllocator());
                            }
                            reference->value.RemoveMember("feed by segment");
                            reference->value.AddMember("feed by segment", feed_by_segment.Move(), instruction.GetAllocator());
                        }
                    }
                }

                reference = value.FindMember("codec");
                if(reference != value.MemberEnd()) {
                    if(!reference->value.IsNull()) {
                        if(reference->value.IsObject()) {
                            for(auto& record : reference->value.GetObject()) {
                                list< URL > feed_url_array(decode_value_by_key< list< URL > >("output", record.value));
                                Value feed_by_segment(kArrayType);
                                for(const auto& url : feed_url_array) {
                                    const Value& proxy(feed_ontology_by_url[url]);
                                    feed_by_segment.PushBack(Value(proxy, instruction.GetAllocator()).Move(), instruction.GetAllocator());
                                }
                                record.value.RemoveMember("feed by segment");
                                record.value.AddMember("feed by segment", feed_by_segment.Move(), instruction.GetAllocator());
                            }
                        }
                    }
                }

                Value feed_array(kArrayType);
                for(auto& record : feed_ontology_by_url) {
                    feed_array.PushBack(record.second.Move(), instruction.GetAllocator());
                }
                instruction.RemoveMember("output feed");
                instruction.AddMember("output feed", feed_array.Move(), instruction.GetAllocator());
            }
            cross_validate_io();
        }
    }
};
void Demultiplex::cross_validate_io() {
    list< URL > input_array(decode_value_by_key< list< URL > >("input", instruction));

    set< URL > input;
    for(auto& url : input_array) {
        input.emplace(url);
    }

    Value::MemberIterator reference = instruction.FindMember("multiplex");
    if(reference != instruction.MemberEnd()) {
        if(reference->value.IsObject()) {
            reference = reference->value.FindMember("codec");
            if(reference != instruction.MemberEnd()) {
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
void Demultiplex::validate_instruction() {
    Action::validate_instruction();

    uint8_t input_phred_offset;
    if(decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, instruction)) {
        if(input_phred_offset > MAX_PHRED_VALUE || input_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("input phred offset out of range " + to_string(input_phred_offset));
        }
    }

    uint8_t output_phred_offset;
    if(decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, instruction)) {
        if(output_phred_offset > MAX_PHRED_VALUE || output_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("output phred offset out of range " + to_string(output_phred_offset));
        }
    }

    validate_codec_group_sanity("multiplex");
    validate_codec_group_sanity("molecular");
    validate_codec_group_sanity("splitseq");
};
void Demultiplex::validate_codec_sanity(Value& value) {
    if(!value.IsNull()) {
        CodecMetric metric(value);
        if(!metric.empty()) {
            metric.apply_barcode_tolerance(value, instruction);
        }
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
};
void Demultiplex::validate_codec_group_sanity(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                validate_codec_sanity(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& decoder : reference->value.GetArray()) {
                    if(!decoder.IsNull()) {
                        validate_codec_sanity(decoder);
                    }
                }
            }
        }
    }
};
void Demultiplex::clean_instruction() {
    instruction.RemoveMember("decoder");
    Action::clean_instruction();
};
void Demultiplex::compile_PG_element() {
    Value PG(kObjectType);

    string buffer;
    if(decode_value_by_key< string >("application name", buffer, instruction)) {
        encode_key_value("ID", buffer, PG, instruction);
    }
    if(decode_value_by_key< string >("application name", buffer, instruction)) {
        encode_key_value("PN", buffer, PG, instruction);
    }
    if(decode_value_by_key< string >("full command", buffer, instruction)) {
        encode_key_value("CL", buffer, PG, instruction);
    }
    if(decode_value_by_key< string >("previous application", buffer, instruction)) {
        encode_key_value("PP", buffer, PG, instruction);
    }
    if(decode_value_by_key< string >("application description", buffer, instruction)) {
        encode_key_value("DS", buffer, PG, instruction);
    }
    if(decode_value_by_key< string >("application version", buffer, instruction)) {
        encode_key_value("VN", buffer, PG, instruction);
    }
    instruction.AddMember(Value("program", instruction.GetAllocator()).Move(), PG.Move(), instruction.GetAllocator());
};

void Interface::load_sub_action(const Value& ontology) {
    string key;
    if(decode_value_by_key< string >("name", key, ontology)) {
        if(action_by_name.find(key) == action_by_name.end()) {
            Action* action(NULL);
            if(key == "demux") {
                action = new Demultiplex(ontology);
            } else {
                action = new Action(ontology);
            }
            action_by_index.push_back(action);
            action_by_name.emplace(make_pair(key, action));
            layout.max_action_name = max(layout.max_action_name, action->name_length());
        } else { throw ConfigurationError("duplicate action " + key); }
    } else { throw InternalError("undefined action"); }
};
ostream& Interface::print_version_element(ostream& o) const {
    CommandLine::print_version_element(o);

    #ifdef ZLIB_VERSION
        o << "zlib " << ZLIB_VERSION << endl;
    #endif

    #ifdef BZIP2_VERSION
        o << "bzlib " << BZIP2_VERSION << endl;
    #endif

    #ifdef XZ_VERSION
        o << "xzlib " << XZ_VERSION << endl;
    #endif

    #ifdef LIBDEFLATE_VERSION
        o << "libdeflate " << LIBDEFLATE_VERSION << endl;
    #endif

    #ifdef RAPIDJSON_VERSION
        o << "rapidjson " << RAPIDJSON_VERSION << endl;
    #endif

    #ifdef HTSLIB_VERSION
        o << "htslib " << HTSLIB_VERSION << endl;
    #endif

    return o;
};
