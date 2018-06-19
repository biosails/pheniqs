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

template<> CodecDistanceMetric* decode_value_by_key(const Value::Ch* key, const Value& container) {
    CodecDistanceMetric* value(NULL);
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                value = new CodecDistanceMetric(reference->value);
            } else { throw ConfigurationError("codec element must be a dictionary"); }
        }
    }
    return value;
};
template<> bool decode_value_by_key< list< CodecDistanceMetric > >(const Value::Ch* key, list< CodecDistanceMetric >& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsArray()) {
            if(!element->value.Empty()) {
                for(const auto& e : element->value.GetArray()) {
                    value.emplace_back(e);
                }
                return true;
            }
        } else { throw ConfigurationError(string(key) + " element must be an array"); }
    }
    return false;
};

void Demultiplex::apply_instruction_manipulation() {
    apply_codec_inheritence();

    embed_codec("multiplex");
    project_codec_group("multiplex");
    infer_codec_group("multiplex");
    enumerate_codec_group("multiplex");
    normalize_codec_group_concentration("multiplex");
    expand_codec_group_url_variable("multiplex");
    expand_codec_group_url_base("multiplex");
    load_input_feed();
    cross_validate_codec_group_io("multiplex");
    load_codec_group_transformation("multiplex");
    load_output_transformation("multiplex");
    embed_codec("molecular");
    project_codec_group("molecular");
    enumerate_codec_group("molecular");
    normalize_codec_group_concentration("molecular");
    load_codec_group_transformation("molecular");

    embed_codec("splitseq");
    project_codec_group("splitseq");
    enumerate_codec_group("splitseq");
    normalize_codec_group_concentration("splitseq");
    load_codec_group_transformation("splitseq");
};
void Demultiplex::load_input_feed() {
    Platform platform(Platform::UNKNOWN);
    decode_value_by_key< Platform >("platform", platform, instruction);

    int32_t buffer_capacity(numeric_limits< int32_t >::max());
    decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, instruction);

    uint8_t input_phred_offset(numeric_limits< uint8_t >::max());
    decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, instruction);

    list< URL > feed_url_array(decode_value_by_key< list< URL > >("input", instruction));
    int32_t input_segment_cardinality(static_cast< int32_t>(feed_url_array.size()));
    encode_key_value("input segment cardinality", input_segment_cardinality, instruction, instruction);

    /*  validate leading_segment_index */
    int32_t leading_segment_index(numeric_limits< int32_t >::max());
    if(decode_value_by_key< int32_t >("leading segment index", leading_segment_index, instruction)) {
        if(leading_segment_index >= input_segment_cardinality) {
            throw ConfigurationError("invalid leading segment index " + to_string(leading_segment_index));
        }
    }

    map< URL, int > feed_resolution;
    for(const auto& url : feed_url_array) {
        ++(feed_resolution[url]);
    }

    unordered_map< URL, Value > feed_ontology_by_url;
    int32_t feed_index(0);
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
    instruction.AddMember("input feed by segment", feed_by_segment.Move(), instruction.GetAllocator());

    Value feed_array(kArrayType);
    for(auto& record : feed_ontology_by_url) {
        feed_array.PushBack(record.second.Move(), instruction.GetAllocator());
    }
    instruction.AddMember("input feed", feed_array.Move(), instruction.GetAllocator());
};
int32_t Demultiplex::inheritence_depth(const string& key, const unordered_map< string, Value* >& node_by_key, Document& document) {
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
    } else { throw ConfigurationError("referencing an undefined base codec " + key); }
    return depth;
};
void Demultiplex::apply_codec_inheritence() {
    unordered_map< string, Value* > codec_by_key;
    Value::MemberIterator reference = instruction.FindMember("codec");
    if(reference != instruction.MemberEnd()) {
        Value& original = reference->value;
        if(!original.IsNull()) {
            if(original.IsObject()) {

                /*  map a copy of each codec by key */
                for(auto& codec_record : original.GetObject()) {
                    if(!codec_record.value.IsNull()) {
                        string codec_key(codec_record.name.GetString(), codec_record.name.GetStringLength());
                        codec_record.value.RemoveMember("depth");
                        codec_by_key.emplace(make_pair(codec_key, &codec_record.value));
                    }
                }

                /* compute the inheritence depth of each codec and keep track of the max depth */
                int32_t max_depth(0);
                for(auto& record : codec_by_key) {
                    max_depth = max(max_depth, inheritence_depth(record.first, codec_by_key, instruction));
                }

                /* apply codec inheritence down the tree */
                int32_t depth(0);
                for(int32_t i(1); i <= max_depth; ++i) {
                    for(auto& record : codec_by_key) {
                        Value* value = record.second;
                        if(decode_value_by_key("depth", depth, *value) && depth == i) {
                            string base_key;
                            if(decode_value_by_key< string >("base", base_key, *value)) {
                                auto base_record = codec_by_key.find(base_key);
                                if(base_record != codec_by_key.end()) {
                                    Value* base = base_record->second;
                                    merge_json_value(*base, *value, instruction);
                                }
                            }
                        }
                    }
                }
            } else { throw ConfigurationError("codec element must be a dictionary"); }
        }
    }
};
void Demultiplex::embed_codec(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember("codec");
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                if(!reference->value.ObjectEmpty()) {

                    /* create a map of globally available codecs */
                    unordered_map< string, Value* > codec_by_key;
                    for(auto& codec_record : reference->value.GetObject()) {
                        if(!codec_record.value.IsNull()) {
                            string key(codec_record.name.GetString(), codec_record.name.GetStringLength());
                            codec_by_key.emplace(make_pair(key, &codec_record.value));
                        }
                    }

                    reference = instruction.FindMember(key);
                    if(reference != instruction.MemberEnd()) {
                        if(!reference->value.IsNull()) {
                            if(reference->value.IsObject()) {
                                string base_key;
                                if(decode_value_by_key< string >("base", base_key, reference->value)) {
                                    auto base_record = codec_by_key.find(base_key);
                                    if(base_record != codec_by_key.end()) {
                                        merge_json_value(*base_record->second, reference->value, instruction);
                                    }
                                }
                                encode_key_value("index", 0, reference->value, instruction);

                            } else if(reference->value.IsArray()) {
                                int32_t index(0);
                                for(auto& element : reference->value.GetArray()) {
                                    if(element.IsObject()) {
                                        string base_key;
                                        if(decode_value_by_key< string >("base", base_key, element)) {
                                            auto base_record = codec_by_key.find(base_key);
                                            if(base_record != codec_by_key.end()) {
                                                merge_json_value(*base_record->second, element, instruction);
                                            }
                                        }
                                        encode_key_value("index", index, element, instruction);
                                        ++index;
                                    } else { throw ConfigurationError("codec element must be a dictionary"); }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};
void Demultiplex::project_codec(Value& value, const Value& default_instruction_codec, const Value& default_instruction_barcode) {
    if(!value.IsNull()) {
        merge_json_value(default_instruction_codec, value, instruction);

        Value default_barcode(kObjectType);
        project_json_value(default_instruction_barcode, value, default_barcode, instruction);

        Value::MemberIterator reference = value.FindMember("barcode");
        if(reference != value.MemberEnd()){
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    for(auto& record : reference->value.GetObject()) {
                        merge_json_value(default_barcode, record.value, instruction);
                    }
                } else { throw ConfigurationError("barcode element must be a dictionary"); }
            }
        }

        reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()){
            merge_json_value(default_barcode, reference->value, instruction);
        } else {
            value.AddMember (
                Value("undetermined", instruction.GetAllocator()).Move(),
                Value(default_barcode, instruction.GetAllocator()).Move(),
                instruction.GetAllocator()
            );
        }
        default_barcode.SetNull();
    }
};
void Demultiplex::project_codec_group(const Value::Ch* key) {
    /* get the base barcode and codec node from the configuration */
    Value default_configuration_codec(kObjectType);
    Value default_configuration_barcode(kObjectType);
    Value::ConstMemberIterator const_reference = ontology.FindMember("projection");
    if(const_reference != ontology.MemberEnd()) {
        const Value& base = const_reference->value;
        if(!base.IsNull()) {
            const_reference = base.FindMember("codec");
            if(const_reference != base.MemberEnd() && !const_reference->value.IsNull()) {
                default_configuration_codec.CopyFrom(const_reference->value, instruction.GetAllocator());
            }
            const_reference = base.FindMember("barcode");
            if(const_reference != base.MemberEnd() && !const_reference->value.IsNull()) {
                default_configuration_barcode.CopyFrom(const_reference->value, instruction.GetAllocator());
            }
        }
    }

    /* project codec attributes from the document root */
    Value default_instruction_codec(kObjectType);
    project_json_value(default_configuration_codec, instruction, default_instruction_codec, instruction);
    default_configuration_codec.SetNull();

    /* project barcode attributes from the instruction root */
    Value default_instruction_barcode(kObjectType);
    project_json_value(default_configuration_barcode, instruction, default_instruction_barcode, instruction);
    default_configuration_barcode.SetNull();

    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                project_codec(reference->value, default_instruction_codec, default_instruction_barcode);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        project_codec(codec, default_instruction_codec, default_instruction_barcode);
                    }
                }
            }
        }
    }
    default_instruction_codec.SetNull();
    default_instruction_barcode.SetNull();
};
bool Demultiplex::infer_PU(const Value::Ch* key, string& value, Value& container, const bool& undetermined) {
    if(decode_value_by_key< string >(key, value, container)) {
        return true;
    } else {
        string PU;
        if(undetermined) {
            PU.assign("undetermined");
        } else {
            list< string > barcode;
            if(decode_value_by_key< list< string > >("segment", barcode, container)) {
                for(auto& segment : barcode) {
                    PU.append(segment);
                }
            }
        }
        if(!PU.empty()) {
            string flowcell_id;
            if(decode_value_by_key< string >("flowcell id", flowcell_id, container)) {
                value.assign(flowcell_id);
                value.push_back(':');

                int32_t flowcell_lane_number;
                if(decode_value_by_key< int32_t >("flowcell lane number", flowcell_lane_number, container)) {
                    value.append(to_string(flowcell_lane_number));
                    value.push_back(':');
                }
            }
            value.append(PU);
            encode_key_value(key, value, container, instruction);
            return true;
        }
    }
    return false;
};
bool Demultiplex::infer_ID(const Value::Ch* key, string& value, Value& container, const bool& undetermined) {
    if(decode_value_by_key< string >(key, value, container)) {
        return true;
    } else if(infer_PU("PU", value, container, undetermined)) {
        encode_key_value(key, value, container, instruction);
        return true;
    } else { return false; }
};
void Demultiplex::enumerate_codec(Value& value) {
    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                int32_t index(1);
                for(auto& record : reference->value.GetObject()) {
                    if(!record.value.IsNull()) {
                        encode_key_value("index", index, record.value, instruction);
                        ++index;
                    }
                }
            }
        }
    }

    reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            encode_key_value("index", 0, reference->value, instruction);
        }
    }
};
void Demultiplex::enumerate_codec_group(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                enumerate_codec(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        enumerate_codec(codec);
                    }
                }
            }
        }
    }
};
void Demultiplex::infer_codec(Value& value) {
    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            set< string > unique_word_id;
            for(auto& record : reference->value.GetObject()) {
                if(!record.value.IsNull()) {
                    string key;
                    if(infer_ID("ID", key, record.value)) {
                        if(!unique_word_id.count(key)) {
                            unique_word_id.emplace(key);
                        } else {
                            string duplicate(record.name.GetString(), record.value.GetStringLength());
                            throw ConfigurationError("duplicate " + duplicate + " barcode");
                        }
                    }
                }
            }
        }
    }
    reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            string key;
            infer_ID("ID", key, reference->value, true);
        }
    }
};
void Demultiplex::infer_codec_group(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                infer_codec(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        infer_codec(codec);
                    }
                }
            }
        }
    }
};
void Demultiplex::normalize_codec_concentration(Value& value) {
    double noise(decode_value_by_key< double >("noise", value));

    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        Value& barcode_dictionary = reference->value;
        if(!barcode_dictionary.IsNull()) {
            if(barcode_dictionary.IsObject()) {
                double total(0);
                for(auto& barcode_record : barcode_dictionary.GetObject()) {
                    double concentration(decode_value_by_key< double >("concentration", barcode_record.value));
                    if(concentration >= 0) {
                        total += concentration;
                    } else { throw ConfigurationError("barcode word concentration must be a positive number");  }
                }

                const double factor((1.0 - noise) / total);
                for(auto& barcode_record : barcode_dictionary.GetObject()) {
                    double concentration(decode_value_by_key< double >("concentration", barcode_record.value));
                    encode_key_value("concentration", concentration * factor, barcode_record.value, instruction);
                }
            }
        }
    }

    reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                encode_key_value("concentration", noise, reference->value, instruction);
            }
        }
    }
};
void Demultiplex::normalize_codec_group_concentration(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                normalize_codec_concentration(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        normalize_codec_concentration(codec);
                    }
                }
            }
        }
    }
};
void Demultiplex::expand_codec_url_variable(Value& value) {
    /*  expand shell variables and standard streams in URL objects */
    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        Value& barcode_dictionary = reference->value;
        if(!barcode_dictionary.IsNull()) {
            if(barcode_dictionary.IsObject()) {
                for(auto& barcode_record : barcode_dictionary.GetObject()) {
                    reference = barcode_record.value.FindMember("output");
                    if(reference != barcode_record.value.MemberEnd()) {
                        Value& original_output_array = reference->value;
                        if(!original_output_array.IsNull() && !original_output_array.Empty()) {
                            if(original_output_array.IsArray()) {
                                Value expanded;
                                Value output_array(kArrayType);
                                for(const auto& url : original_output_array.GetArray()) {
                                    if(expand_file_url(url, expanded, instruction, IoDirection::OUT)) {
                                        output_array.PushBack(expanded.Move(), instruction.GetAllocator());
                                    }
                                }
                                barcode_record.value.RemoveMember("output");
                                barcode_record.value.AddMember(Value("output", instruction.GetAllocator()).Move(), output_array.Move(), instruction.GetAllocator());
                            }
                        }
                    }
                }
            }
        }
    }
};
void Demultiplex::expand_codec_group_url_variable(const Value::Ch* key) {
    Value expanded;
    Value::MemberIterator reference = instruction.FindMember("base input url");
    if(reference != instruction.MemberEnd() && !reference->value.IsNull()) {
        if(expand_directory_url(reference->value, expanded, instruction)) {
            instruction.RemoveMember("base input url");
            instruction.AddMember(Value("base input url", instruction.GetAllocator()).Move(), expanded.Move(), instruction.GetAllocator());
        }
    }

    reference = instruction.FindMember("base output url");
    if(reference != instruction.MemberEnd() && !reference->value.IsNull()) {
        if(expand_directory_url(reference->value, expanded, instruction)) {
            instruction.RemoveMember("base output url");
            instruction.AddMember(Value("base output url", instruction.GetAllocator()).Move(), expanded.Move(), instruction.GetAllocator());
        }
    }

    reference = instruction.FindMember("input url");
    if(reference != instruction.MemberEnd() && !reference->value.IsNull()) {
        if(expand_file_url(reference->value, expanded, instruction, IoDirection::IN)) {
            instruction.RemoveMember("input url");
            instruction.AddMember(Value("input url", instruction.GetAllocator()).Move(), expanded.Move(), instruction.GetAllocator());
        }
    }

    reference = instruction.FindMember("input");
    if(reference != instruction.MemberEnd() && !reference->value.IsNull()) {
        if(reference->value.IsArray() && !reference->value.Empty()) {
            Value array(kArrayType);
            for(const auto& url : reference->value.GetArray()) {
                if(expand_file_url(url, expanded, instruction, IoDirection::IN)){
                    array.PushBack(expanded.Move(), instruction.GetAllocator());
                }
            }
            instruction.RemoveMember("input");
            instruction.AddMember(Value("input", instruction.GetAllocator()).Move(), array.Move(), instruction.GetAllocator());
        }
    }

    reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                expand_codec_url_variable(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        expand_codec_url_variable(codec);
                    }
                }
            }
        }
    }
};
void Demultiplex::expand_codec_url_base(Value& value) {
    URL base_output_url;
    if(decode_directory_url_by_key("base output url", base_output_url, instruction)) {
        Value::MemberIterator reference = value.FindMember("barcode");
        if(reference != value.MemberEnd()) {
            Value& barcode_dictionary = reference->value;
            if(!barcode_dictionary.IsNull()) {
                if(barcode_dictionary.IsObject()) {
                    for(auto& barcode_record : barcode_dictionary.GetObject()) {
                        list< URL > output;
                        if(decode_value_by_key< list< URL > >("output", output, barcode_record.value)) {
                            bool output_disabled(true);
                            for(auto& url : output) {
                                url.relocate(base_output_url);
                                output_disabled = output_disabled && url.is_null();
                            }
                            encode_key_value("output", output, barcode_record.value, instruction);
                            encode_key_value("output disabled", output_disabled, barcode_record.value, instruction);
                        }
                    }
                }
            }
        }

        reference = value.FindMember("undetermined");
        if(reference != value.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    list< URL > output;
                    if(decode_value_by_key< list< URL > >("output", output, reference->value)) {
                        bool output_disabled(true);
                        for(auto& url : output) {
                            url.relocate(base_output_url);
                            output_disabled = output_disabled && url.is_null();
                        }
                        encode_key_value("output", output, reference->value, instruction);
                        encode_key_value("output disabled", output_disabled, reference->value, instruction);
                    }
                }
            }
        }
    }
};
void Demultiplex::expand_codec_group_url_base(const Value::Ch* key) {
    URL base_input_url;
    if(decode_directory_url_by_key("base input url", base_input_url, instruction)) {
        list< URL > input;
        if(decode_value_by_key< list< URL > >("input", input, instruction)) {
            for(auto& url : input) {
                url.relocate(base_input_url);
            }
            encode_key_value("input", input, instruction, instruction);
        }

        URL input_url;
        if(decode_value_by_key< URL >("input url", input_url, instruction)) {
            input_url.relocate(base_input_url);
            encode_key_value("input url", input_url, instruction, instruction);
        }
    }

    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                expand_codec_url_base(reference->value);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        expand_codec_url_base(codec);
                    }
                }
            }
        }
    }

};
void Demultiplex::cross_validate_codec_io(Value& value, const set< URL >& input) {
    /*  verify no URL is used for both input and output */
    Value::MemberIterator reference = value.FindMember("output feed");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull() && !reference->value.Empty()) {
            for(auto& element : reference->value.GetArray()) {
                URL url;
                if(decode_value_by_key< URL >("url", url, element)) {
                    if(input.count(url) > 0) {
                        throw ConfigurationError("URL " + string(url) + " is used for both input and output");
                    }
                }
            }
        }
    }
};
void Demultiplex::cross_validate_codec_group_io(const Value::Ch* key) {
    Value::MemberIterator reference = instruction.FindMember("input feed");
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull() && !reference->value.Empty()) {
            set< URL > input;
            for(auto& element : reference->value.GetArray()) {
                URL url;
                if(decode_value_by_key< URL >("url", url, element)) {
                    input.emplace(url);
                }
            }

            reference = instruction.FindMember(key);
            if(reference != instruction.MemberEnd()) {
                if(!reference->value.IsNull()) {
                    if(reference->value.IsObject()) {
                        cross_validate_codec_io(reference->value, input);
                    } else if(reference->value.IsArray()) {
                        for(auto& codec : reference->value.GetArray()) {
                            if(!codec.IsNull()) {
                                cross_validate_codec_io(codec, input);
                            }
                        }
                    }
                }
            }
        }
    }
};
void Demultiplex::pad_codec_output_url(Value& value, const int32_t& output_segment_cardinality) {
    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        Value& barcode_dictionary = reference->value;
        if(!barcode_dictionary.IsNull()) {
            if(barcode_dictionary.IsObject()) {
                for(auto& barcode_record : barcode_dictionary.GetObject()) {
                    list< URL > array;
                    if(decode_value_by_key< list< URL > >("output", array, barcode_record.value)) {
                        if(!array.empty()) {
                            if(static_cast< int32_t >(array.size()) != output_segment_cardinality) {
                                if(array.size() == 1) {
                                    while(static_cast< int32_t >(array.size()) < output_segment_cardinality) {
                                        array.push_back(array.front());
                                    }
                                    encode_key_value("output", array, barcode_record.value, instruction);
                                } else {
                                    string barcode_key(barcode_record.name.GetString(), barcode_record.name.GetStringLength());
                                    throw ConfigurationError("incorrect number of output URLs in channel " + barcode_key);
                                }
                            }
                        }
                    }
                    encode_key_value("TC", output_segment_cardinality, barcode_record.value, instruction);
                }
            }
        }
    }

    reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                Value& undetermined(reference->value);
                list< URL > array;
                if(decode_value_by_key< list< URL > >("output", array, undetermined)) {
                    if(!array.empty()) {
                        if(static_cast< int32_t >(array.size()) != output_segment_cardinality) {
                            if(array.size() == 1) {
                                while(static_cast< int32_t >(array.size()) < output_segment_cardinality) {
                                    array.push_back(array.front());
                                }
                                encode_key_value("output", array, undetermined, instruction);
                            } else {
                                throw ConfigurationError("incorrect number of output URLs in undetermined channel");
                            }
                        }
                    }
                }
                encode_key_value("TC", output_segment_cardinality, undetermined, instruction);
            }
        }
    }
};
void Demultiplex::pad_codec_group_output_url(const Value::Ch* key) {
    int32_t segment_cardinality;
    if(decode_value_by_key< int32_t >("output segment cardinality", segment_cardinality, instruction)) {
        Value::MemberIterator reference = instruction.FindMember(key);
        if(reference != instruction.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    pad_codec_output_url(reference->value, segment_cardinality);
                } else if(reference->value.IsArray()) {
                    for(auto& codec : reference->value.GetArray()) {
                        if(!codec.IsNull()) {
                            pad_codec_output_url(codec, segment_cardinality);
                        }
                    }
                }
            }
        }
    }
};
void Demultiplex::load_codec_output_feed(Value& value, const Platform& platform, const int32_t& buffer_capacity, const uint8_t& phred_offset) {
    unordered_map< URL, unordered_map< int32_t, int > > feed_resolution;

    Value::MemberIterator reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            int32_t index(decode_value_by_key< int32_t >("index", reference->value));
            list< URL > feed_url_array(decode_value_by_key< list< URL > >("output", reference->value));
            for(auto& url : feed_url_array) {
                ++(feed_resolution[url][index]);
            }
        }
    }

    reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                for(auto& record : reference->value.GetObject()) {
                    int32_t index(decode_value_by_key< int32_t >("index", record.value));
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
                    reference->value.AddMember("feed by segment", feed_by_segment.Move(), instruction.GetAllocator());
                }
            }
        }

        reference = value.FindMember("barcode");
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
                        record.value.AddMember("feed by segment", feed_by_segment.Move(), instruction.GetAllocator());
                    }
                }
            }
        }

        Value feed_array(kArrayType);
        for(auto& record : feed_ontology_by_url) {
            feed_array.PushBack(record.second.Move(), instruction.GetAllocator());
        }
        instruction.AddMember("output feed", feed_array.Move(), instruction.GetAllocator());
    }
};
void Demultiplex::load_codec_group_output_feed(const Value::Ch* key) {
    Platform platform(Platform::UNKNOWN);
    decode_value_by_key< Platform >("platform", platform, instruction);

    int32_t buffer_capacity(numeric_limits< int32_t >::max());
    decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, instruction);

    uint8_t output_phred_offset(numeric_limits< uint8_t >::max());
    decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, instruction);

    Value::MemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                load_codec_output_feed(reference->value, platform, buffer_capacity, output_phred_offset);
            } else if(reference->value.IsArray()) {
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        load_codec_output_feed(codec, platform, buffer_capacity, output_phred_offset);
                    }
                }
            }
        }
    }
};
void Demultiplex::load_output_transformation(const Value::Ch* key) {
    int32_t input_segment_cardinality(numeric_limits< int32_t >::max());
    if(decode_value_by_key< int32_t >("input segment cardinality", input_segment_cardinality, instruction)) {
        Rule rule(decode_value_by_key< Rule >("template", instruction));
        for(auto& token : rule.token_array) {
            if(!(token.input_segment_index < input_segment_cardinality)) {
                throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
            }
        }
        encode_key_value("output segment cardinality", rule.output_segment_cardinality, instruction, instruction);
        pad_codec_group_output_url(key);
        load_codec_group_output_feed(key);
    }
};
void Demultiplex::load_codec_transformation(Value& value, const int32_t& input_segment_cardinality) {
    Rule rule(decode_value_by_key< Rule >("template", value));
    for(auto& token : rule.token_array) {
        if(!(token.input_segment_index < input_segment_cardinality)) {
            throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
        }
    }

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

    Value::MemberIterator reference = value.FindMember("barcode");
    if(reference != value.MemberEnd()) {
        Value& barcode_dictionary = reference->value;
        if(!barcode_dictionary.IsNull()) {
            if(barcode_dictionary.IsObject()) {
                for(auto& barcode_record : barcode_dictionary.GetObject()) {
                    encode_key_value("segment cardinality", rule.output_segment_cardinality, barcode_record.value, instruction);
                }
            }
        }
    }

    reference = value.FindMember("undetermined");
    if(reference != value.MemberEnd()) {
        if(!reference->value.IsNull()) {
            Value& undetermined(reference->value);

            Value segment(kArrayType);
            for(size_t i = 0; i < barcode_length.size(); ++i) {
                string sequence(barcode_length[i], '=');
                segment.PushBack(Value(sequence.c_str(), sequence.size(), instruction.GetAllocator()).Move(), instruction.GetAllocator());
            }
            undetermined.RemoveMember("segment");
            undetermined.AddMember(Value("segment", instruction.GetAllocator()).Move(), segment.Move(), instruction.GetAllocator());

            encode_key_value("segment cardinality", rule.output_segment_cardinality, undetermined, instruction);
        }
    }
};
void Demultiplex::load_codec_group_transformation(const Value::Ch* key) {
    int32_t input_segment_cardinality;
    if(decode_value_by_key< int32_t >("input segment cardinality", input_segment_cardinality, instruction)) {
        Value::MemberIterator reference = instruction.FindMember(key);
        if(reference != instruction.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsObject()) {
                    load_codec_transformation(reference->value, input_segment_cardinality);
                } else if(reference->value.IsArray()) {
                    for(auto& codec : reference->value.GetArray()) {
                        if(!codec.IsNull()) {
                            load_codec_transformation(codec, input_segment_cardinality);
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
        CodecDistanceMetric metric(value);
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
                for(auto& codec : reference->value.GetArray()) {
                    if(!codec.IsNull()) {
                        validate_codec_sanity(codec);
                    }
                }
            }
        }
    }
};
void Demultiplex::clean_instruction() {
    instruction.RemoveMember("codec");
    Action::clean_instruction();
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
