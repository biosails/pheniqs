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

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "environment.h"

/*  Environment */
Environment::Environment(const int argc, const char** argv) :
    interface(argc, argv),
    _program_action(interface.get_selected_action()),
    _help_only(interface.help_triggered()),
    _version_only(interface.version_triggered()),
    _validate_only(false),
    _lint_only(false),
    _display_distance(false),
    _instruction(kObjectType) {

    if(!_help_only && !_version_only) {
        switch (_program_action) {
            case ProgramAction::DEMULTIPLEX: {
                load_interface_instruction();
                decode_value_by_key< bool >("validate only", _validate_only, _instruction);
                decode_value_by_key< bool >("lint only", _lint_only, _instruction);
                decode_value_by_key< bool >("display distance", _display_distance, _instruction);
                break;
            };

            case ProgramAction::QUALITY: {
                break;
            };

            default: break;
        }
    }
};
void Environment::load_interface_instruction() {
    const Document& original = interface.instruction();
    Value::ConstMemberIterator collection;
    Document base;

    /*  extract a base read group from the root element */
    transcode_head_RG_atom(original, base, base);

    /*  overlay any explicitly defined read groups on top of the default */
    map< string, Value* > read_group_node_by_id;
    collection = original.FindMember("read group");
    if(collection != original.MemberEnd()) {
        if(collection->value.IsArray() && !collection->value.Empty()) {
            for(auto& element : collection->value.GetArray()) {
                string key;
                decode_value_by_key< string >("ID", key, element);
                if(!key.empty()) {
                    Value* implicit_read_group = new Value();
                    auto record = read_group_node_by_id.find(key);
                    if(record == read_group_node_by_id.end()) {
                        merge_json_value(base, element, *implicit_read_group, _instruction);
                        read_group_node_by_id.emplace(make_pair(key, implicit_read_group));
                    } else {
                        merge_json_value(*record->second, element, *implicit_read_group, _instruction);
                        delete record->second;
                        read_group_node_by_id[key] = implicit_read_group;
                    }
                } else { throw ConfigurationError("read group is missing an ID"); }
            }
        } else { throw ConfigurationError("read group element must be an array"); }
    }

    /*  extract a base channel from the root element */
    transcode_channel_specification(original, base, base);

    /*  if a matching read group has been defined overlay the read group 
        on top of the default channel and than the channel on top of that.
        otherwise overlay every channel on top of the default 

        also make sure there are no two channels with the same read group ID */
    collection = original.FindMember("channel");
    if(collection != original.MemberEnd()) {
        set< string > unique_channel_id;
        if(collection->value.IsArray() && !collection->value.Empty()) {
            uint64_t channel_index(0);
            uint64_t undetermined_count(0);
            bool undetermined(false);
            Value channel_array(kArrayType);
            for(auto& element : collection->value.GetArray()) {
                string key;
                decode_value_by_key< string >("RG", key, element);
                if(!key.empty()) {
                    if(!unique_channel_id.count(key)) {
                        unique_channel_id.emplace(key);

                        Value implicit_channel;
                        auto record = read_group_node_by_id.find(key);
                        if(record == read_group_node_by_id.end()) {
                            merge_json_value(base, element, implicit_channel, _instruction);
                        } else {
                            Value annotated_channel;
                            merge_json_value(base, *record->second, annotated_channel, _instruction);
                            merge_json_value(annotated_channel, element, implicit_channel, _instruction);
                            annotated_channel.SetNull();
                        }

                        /* remove any barcode and concentration declared on the undetermined channel */
                        if(decode_value_by_key< bool >("undetermined", undetermined, implicit_channel)) {
                            if(undetermined) {
                                implicit_channel.RemoveMember("barcode");
                                implicit_channel.RemoveMember("concentration");
                                undetermined_count++;
                            }
                        }

                        encode_key_value("index", channel_index++, implicit_channel, _instruction);
                        channel_array.PushBack(implicit_channel.Move(), _instruction.GetAllocator());
                    } else { throw ConfigurationError("channel " + key + " at position " + to_string(channel_index) + " is already defined"); }
                } else { throw ConfigurationError("channel is missing a read group ID"); }
            }

            if(undetermined_count > 1) {
                throw ConfigurationError("only one undetermined channel can be declared");

            } else if(channel_index == 1 && undetermined_count == 0) {
                /*  if only one channel is defined and it has no barcode defined, mark it as undetermined */
                for(auto& element : channel_array.GetArray()) {
                    Barcode barcode;
                    if(!decode_value_by_key< Barcode >("barcode", barcode, element)) {
                        encode_key_value("undetermined", true, element, _instruction);
                    }
                }
            }
            _instruction.AddMember("channel", channel_array.Move(), _instruction.GetAllocator());
        } else { throw ConfigurationError("channel element must be an array"); }
    }
    base.SetNull();

    /*  encode the read group array */
    if(!read_group_node_by_id.empty()) {
        Value read_group_array(kArrayType);
        for(const auto& record : read_group_node_by_id) {
            read_group_array.PushBack(record.second->Move(), _instruction.GetAllocator());
        }
        _instruction.AddMember("read group", read_group_array.Move(), _instruction.GetAllocator());
    }

    /*  copy everything except the read group and channel arrays */
    for(auto& record : original.GetObject()) {
        string key(record.name.GetString(), record.name.GetStringLength());
        if(key != "read group" && key != "channel") {
            Value value;
            value.CopyFrom(record.value, _instruction.GetAllocator());
            Value name(record.name, _instruction.GetAllocator());
            _instruction.AddMember(name.Move(), value.Move(), _instruction.GetAllocator());
        }
    }
    validate_global_parameters();
    apply_url_base();
    load_concentration_prior();
    load_input_feed_array();
    load_transformation_array();
    cross_validate_io();
    load_pheniqs_pg();
};
void Environment::validate_global_parameters() {
    uint8_t input_phred_offset(numeric_limits< uint8_t >::max());
    if(decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, _instruction)) {
        if(input_phred_offset > MAX_PHRED_VALUE || input_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("input phred offset out of range " + to_string(input_phred_offset));
        }
    }

    uint8_t output_phred_offset(numeric_limits< uint8_t >::max());
    if(decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, _instruction)) {
        if(output_phred_offset > MAX_PHRED_VALUE || output_phred_offset < MIN_PHRED_VALUE) {
            throw ConfigurationError("output phred offset out of range " + to_string(output_phred_offset));
        }
    }

    switch (_program_action) {
        case ProgramAction::DEMULTIPLEX: {
            double multiplex_confidence(numeric_limits< double >::infinity());
            if(decode_value_by_key< double >("multiplex confidence", multiplex_confidence, _instruction)) {
                if(multiplex_confidence < 0 || multiplex_confidence > 1) {
                    throw ConfigurationError("multiplex confidence value " + to_string(multiplex_confidence) + " not between 0 and 1");
                }
            }
            double multiplex_noise(numeric_limits< double >::infinity());
            if(decode_value_by_key< double >("multiplex noise", multiplex_noise, _instruction)) {
                if(multiplex_noise < 0 || multiplex_noise > 1) {
                    throw ConfigurationError("multiplex noise value " + to_string(multiplex_noise) + " not between 0 and 1");
                }
            }
            break;
        };
        case ProgramAction::QUALITY: {
            break;
        };
        default: break;
    }
};
void Environment::apply_url_base() {
    URL base_input_url;
    if(decode_directory_url_by_key("base input url", base_input_url, _instruction)) {
        list< URL > input;
        if(decode_file_url_list_by_key("input", input, _instruction, IoDirection::IN)) {
            for(auto& url : input) {
                url.relocate(base_input_url);
            }
            encode_key_value("input", input, _instruction, _instruction);
        }
    }

    URL base_output_url;
    if(decode_directory_url_by_key("base output url", base_output_url, _instruction)) {
        Value::MemberIterator collection = _instruction.FindMember("channel");
        if(collection != _instruction.MemberEnd()) {
            if(!collection->value.IsNull() && !collection->value.Empty()) {
                for(auto& element : collection->value.GetArray()) {
                    list< URL > output;
                    if(decode_file_url_list_by_key("output", output, element, IoDirection::OUT)) {
                        for(auto& url : output) {
                            url.relocate(base_output_url);
                        }
                        encode_key_value("output", output, element, _instruction);
                    }
                }
            }
        }
    }
};
void Environment::load_concentration_prior() {
    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            double total_concentration(0);
            uint64_t concentration_count(0);
            uint64_t concentration_specified(0);

            for(auto& element : collection->value.GetArray()) {
                bool undetermined(false);
                decode_value_by_key< bool >("undetermined", undetermined, element);
                if(!undetermined) {
                    concentration_count++;
                    double concentration(0);
                    if(decode_value_by_key< double >("concentration", concentration, element)) {
                        if(concentration >= 0) {
                            concentration_specified++;
                            total_concentration += concentration;
                        } else { throw ConfigurationError("channel concentration must be a positive number");  }
                    }
                } else { element.RemoveMember("concentration"); }
            }

            double multiplex_noise(0);
            decode_value_by_key< double >("multiplex noise", multiplex_noise, _instruction);
            if(concentration_specified > 0) {
                if(concentration_specified == concentration_count) {
                    const double factor((1.0 - multiplex_noise) / total_concentration);
                    for(auto& element : collection->value.GetArray()) {
                        bool undetermined(false);
                        decode_value_by_key< bool >("undetermined", undetermined, element);
                        if(!undetermined) {
                            double concentration(0);
                            if(decode_value_by_key< double >("concentration", concentration, element)) {
                                encode_key_value("concentration", concentration * factor, element, _instruction);
                            }
                        }
                    }
                } else { throw ConfigurationError("inconsistent channel concentration specification"); }
            } else {
                /* if no concentrations were given assume a uniform distribution */
                const double factor((1.0 - multiplex_noise) / concentration_count);
                for(auto& element : collection->value.GetArray()) {
                    bool undetermined(false);
                    decode_value_by_key< bool >("undetermined", undetermined, element);
                    if(!undetermined) {
                        encode_key_value("concentration", factor, element, _instruction);
                    }
                }
            }
        }
    }
};
bool Environment::load_input_feed_array() {
    list< URL > input_url_array;
    if(decode_file_url_list_by_key("input", input_url_array, _instruction, IoDirection::IN)) {
        uint64_t input_segment_cardinality = input_url_array.size();
        encode_key_value("input segment cardinality", input_segment_cardinality, _instruction, _instruction);

        /*  validate leading_segment_index */
        uint64_t leading_segment_index(numeric_limits< uint64_t >::max());
        if(decode_value_by_key< uint64_t >("leading segment index", leading_segment_index, _instruction)) {
            if(leading_segment_index >= input_segment_cardinality) {
                throw ConfigurationError("invalid leading segment index " + to_string(leading_segment_index));
            }
        }

        map< URL, int > reference;
        for(const auto& url : input_url_array) {
            reference[url]++;
        }

        Platform platform(Platform::UNKNOWN);
        decode_value_by_key< Platform >("platform", platform, _instruction);

        int32_t buffer_capacity(numeric_limits< int32_t >::max());
        decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, _instruction);

        uint8_t input_phred_offset(numeric_limits< uint8_t >::max());
        decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, _instruction);

        uint64_t feed_index(0);
        Value feed_array(kArrayType);
        for(const auto& record : reference) {
            int resolution(record.second);
            Value specification(kObjectType);
            encode_key_value("index", feed_index++, specification, _instruction);
            encode_key_value("url", record.first, specification, _instruction);
            encode_key_value("direction", IoDirection::IN, specification, _instruction);
            encode_key_value("platform", platform, specification, _instruction);
            encode_key_value("capacity", buffer_capacity * resolution, specification, _instruction);
            encode_key_value("resolution", resolution, specification, _instruction);
            encode_key_value("phred offset", input_phred_offset, specification, _instruction);
            feed_array.PushBack(specification.Move(), _instruction.GetAllocator());
        }
        _instruction.AddMember("input feed", feed_array.Move(), _instruction.GetAllocator());
        return true;
    }
    return false;
};
void Environment::load_transformation_array() {
    uint64_t input_segment_cardinality(numeric_limits< uint64_t >::max());
    if(decode_value_by_key< uint64_t >("input segment cardinality", input_segment_cardinality, _instruction)) {
        vector< Token > token_array;
        if(decode_value_by_key< vector< Token > >("token", token_array, _instruction)) {
            for(auto& token : token_array) {
                if(!(token.input_segment_index < input_segment_cardinality)) {
                    throw ConfigurationError("invalid input feed reference " + to_string(token.input_segment_index) + " in token " + to_string(token.index));
                }
            }

            list< Transform > template_transform_array;
            if(decode_transform_array_by_key("template", template_transform_array, _instruction, token_array)) {
                uint64_t output_segment_cardinality(template_transform_array.size());
                encode_key_value("output segment cardinality", output_segment_cardinality, _instruction, _instruction);
                pad_output_url_array(output_segment_cardinality);
                load_output_feed_array();
            }

            // load multiplex barcode transform array
            list< Transform > multiplex_barcode_transform_array;
            if(decode_transform_array_by_key("multiplex barcode", multiplex_barcode_transform_array, _instruction, token_array)) {
                uint64_t multiplex_segment_cardinality(multiplex_barcode_transform_array.size());
                uint64_t concatenated_multiplex_barcode_length(0);
                vector< uint64_t > multiplex_barcode_length(multiplex_segment_cardinality, 0);

                for(auto& transform : multiplex_barcode_transform_array) {
                    if(transform.token.constant()) {
                        if(!transform.token.empty()) {
                            multiplex_barcode_length[transform.output_segment_index] += transform.token.length();
                            concatenated_multiplex_barcode_length += transform.token.length();
                        } else {
                            throw ConfigurationError("multiplex barcode token " + string(transform.token) + " is empty");
                        }
                    } else {
                        throw ConfigurationError("multiplex barcode token " + string(transform.token) + " is not fixed width");
                    }
                }

                double multiplex_noise(0);
                decode_value_by_key< double >("multiplex noise", multiplex_noise, _instruction);
                double random_multiplex_barcode_probability(1.0 / double(pow(4, (concatenated_multiplex_barcode_length))));
                double adjusted_multiplex_noise_probability(multiplex_noise * random_multiplex_barcode_probability);

                encode_key_value("multiplex segment cardinality", multiplex_segment_cardinality, _instruction, _instruction);
                encode_key_value("concatenated multiplex barcode length", concatenated_multiplex_barcode_length, _instruction, _instruction);
                encode_key_value("multiplex barcode length", multiplex_barcode_length, _instruction, _instruction);
                encode_key_value("random multiplex barcode probability", random_multiplex_barcode_probability, _instruction, _instruction);
                encode_key_value("adjusted multiplex noise probability", adjusted_multiplex_noise_probability, _instruction, _instruction);

                load_undetermined_barcode(multiplex_barcode_length);
                load_multiplex_barcode_distance_metric(multiplex_barcode_length);
            }

            // load molecular barcode transform array
            list< Transform > molecular_barcode_transform_array;
            if(decode_transform_array_by_key("molecular barcode", molecular_barcode_transform_array, _instruction, token_array)) {
                uint64_t molecular_segment_cardinality(molecular_barcode_transform_array.size());
                uint64_t concatenated_molecular_barcode_length(0);
                vector< uint64_t > molecular_barcode_length(molecular_segment_cardinality, 0);

                for(auto& transform : molecular_barcode_transform_array) {
                    if(transform.token.constant()) {
                        if(!transform.token.empty()) {
                            molecular_barcode_length[transform.output_segment_index] += transform.token.length();
                            concatenated_molecular_barcode_length += transform.token.length();
                        } else {
                            throw ConfigurationError("molecular barcode token " + string(transform.token) + " is empty");
                        }
                    } else {
                        throw ConfigurationError("molecular barcode token " + string(transform.token) + " is not fixed width");
                    }
                }

                encode_key_value("molecular segment cardinality", molecular_segment_cardinality, _instruction, _instruction);
                encode_key_value("concatenated molecular barcode length", concatenated_molecular_barcode_length, _instruction, _instruction);
                encode_key_value("molecular barcode length", molecular_barcode_length, _instruction, _instruction);
            }


        }
    }
};
void Environment::pad_output_url_array(const uint64_t& output_segment_cardinality) {
    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            for(auto& element : collection->value.GetArray()) {
                list< URL > array;
                if(decode_file_url_list_by_key("output", array, element, IoDirection::OUT)) {
                    if(!array.empty()) {
                        if(array.size() != output_segment_cardinality) {
                            if(array.size() == 1) {
                                while(array.size() < output_segment_cardinality) {
                                    array.push_back(array.front());
                                }
                            } else { 
                                uint64_t channel_index;
                                decode_value_by_key< uint64_t >("index", channel_index, element);
                                throw ConfigurationError("incorrect number of output URLs in channel " + to_string(channel_index)); 
                            }
                        }
                        encode_key_value("output", array, element, _instruction);
                        encode_key_value("TC", output_segment_cardinality, element, _instruction);
                    }
                }
            }
        }
    }
};
bool Environment::load_output_feed_array() {
    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            unordered_map< URL, unordered_map< uint64_t, int > > reference;
            for(auto& element : collection->value.GetArray()) {
                uint64_t channel_index;
                if(decode_value_by_key< uint64_t >("index", channel_index, element)) {
                    list< URL > output;
                    if(decode_file_url_list_by_key("output", output, element, IoDirection::OUT)) {
                        for(auto& url : output) {
                            reference[url][channel_index]++;
                        }
                    }
                }
            }

            Platform platform(Platform::UNKNOWN);
            decode_value_by_key< Platform >("platform", platform, _instruction);

            int32_t buffer_capacity(numeric_limits< int32_t >::max());
            decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, _instruction);

            uint8_t output_phred_offset(numeric_limits< uint8_t >::max());
            decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, _instruction);

            uint64_t feed_index(0);
            Value feed_array(kArrayType);
            for(const auto& url_record : reference) {
                const URL& url = url_record.first;
                int resolution(0);
                for(const auto& channel_record : url_record.second) {
                    if(resolution == 0) {
                        resolution = channel_record.second;
                    } else if(resolution != channel_record.second) {
                        throw ConfigurationError("inconsistent resolution for " + string(url));
                    }
                }
                Value specification(kObjectType);
                encode_key_value("index", feed_index++, specification, _instruction);
                encode_key_value("url", url, specification, _instruction);
                encode_key_value("direction", IoDirection::OUT, specification, _instruction);
                encode_key_value("platform", platform, specification, _instruction);
                encode_key_value("capacity", buffer_capacity * resolution, specification, _instruction);
                encode_key_value("resolution", resolution, specification, _instruction);
                encode_key_value("phred offset", output_phred_offset, specification, _instruction);
                feed_array.PushBack(specification.Move(), _instruction.GetAllocator());
            }
            _instruction.AddMember("output feed", feed_array.Move(), _instruction.GetAllocator());
            return true;
        }
    }
    return false;
};
void Environment::load_undetermined_barcode(const vector< uint64_t >& multiplex_barcode_length) {
    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            for(auto& element : collection->value.GetArray()) {
                bool undetermined(false);
                if(decode_value_by_key< bool >("undetermined", undetermined, element) && undetermined) {
                    Barcode barcode(multiplex_barcode_length.size());
                    for(size_t i = 0; i < multiplex_barcode_length.size(); i++) {
                        string pattern(multiplex_barcode_length[i], '=');
                        barcode.fill(i, pattern.c_str(), pattern.length());
                    }
                    encode_key_value("barcode", barcode, element, _instruction);
                    break;
                }
            }
        }
    }
};
void Environment::load_multiplex_barcode_distance_metric(const vector< uint64_t >& multiplex_barcode_length) {
    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            _multiplex_barcode_set_distance.resize(multiplex_barcode_length.size());
            for(auto& element : collection->value.GetArray()) {
                bool undetermined(false);
                decode_value_by_key< bool >("undetermined", undetermined, element);
                if(!undetermined) {
                    Barcode multiplex_barcode;
                    if(decode_value_by_key< Barcode >("barcode", multiplex_barcode, element)) {
                        if(!multiplex_barcode.empty()) {
                            for(size_t i = 0; i < multiplex_barcode_length.size(); i++) {
                                if(multiplex_barcode.size(i) == multiplex_barcode_length[i]) {
                                    _multiplex_barcode_set_distance[i].add(multiplex_barcode.iupac_ambiguity(i));
                                } else {
                                    uint64_t channel_index;
                                    decode_value_by_key< uint64_t >("index", channel_index, element);
                                    string message("multiplex barcode ");
                                    message.append(to_string(i));
                                    message.append(" in channel ");
                                    message.append(to_string(channel_index));
                                    message.append(" is ");
                                    message.append(to_string(multiplex_barcode.size(i)));
                                    message.append(" nucleotides long, expected ");
                                    message.append(to_string(multiplex_barcode_length[i]));
                                    throw ConfigurationError(message);
                                }
                            }
                            _multiplex_barcode_distance.add(multiplex_barcode.iupac_ambiguity());
                        }
                    }
                }
            }

            for(auto& metric : _multiplex_barcode_set_distance) {
                metric.load();
            }
            load_barcode_tolerance(_multiplex_barcode_set_distance);

            _multiplex_barcode_distance.load();
        }
    }
};
void Environment::load_barcode_tolerance(const vector< BarcodeDistanceMetric >& multiplex_barcode_set_distance) {
    vector< uint8_t > multiplex_barcode_tolerance;
    if(decode_value_by_key< vector< uint8_t > >("multiplex barcode tolerance", multiplex_barcode_tolerance, _instruction)) {
        if(multiplex_barcode_tolerance.size() == multiplex_barcode_set_distance.size()) {
            for(size_t i = 0; i < multiplex_barcode_set_distance.size(); i++) {
                if(multiplex_barcode_tolerance[i] > multiplex_barcode_set_distance[i].shannon_bound()) {
                    throw ConfigurationError(
                        "multiplex barcode tolerance for segment " + 
                        to_string(i) + 
                        " is higher than the shannon bound " + 
                        to_string(multiplex_barcode_set_distance[i].shannon_bound())
                    );
                }
            }
        } else { throw ConfigurationError("inconsistent multiplex barcode tolerance specification"); }
    } else {
        multiplex_barcode_tolerance.resize(multiplex_barcode_set_distance.size());
        for(size_t i = 0; i < multiplex_barcode_set_distance.size(); i++) {
            multiplex_barcode_tolerance[i] = multiplex_barcode_set_distance[i].shannon_bound();
        }
    }
    encode_key_value("multiplex barcode tolerance", multiplex_barcode_tolerance, _instruction, _instruction);

    Value::MemberIterator collection = _instruction.FindMember("channel");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            for(auto& element : collection->value.GetArray()) {
                encode_key_value("multiplex barcode tolerance", multiplex_barcode_tolerance, element, _instruction);
            }
        }
    }
};
void Environment::cross_validate_io() {
    /*  verify no URL is used for both input and output */
    URL url;
    Value::ConstMemberIterator collection;
    
    collection = _instruction.FindMember("input feed");
    if(collection != _instruction.MemberEnd()) {
        if(!collection->value.IsNull() && !collection->value.Empty()) {
            set< URL > input_feed_url;
            for(auto& element : collection->value.GetArray()) {
                if(decode_file_url_by_key("url", url, IoDirection::IN, element)) {
                    input_feed_url.emplace(url);
                }
            }

            collection = _instruction.FindMember("output feed");
            if(collection != _instruction.MemberEnd()) {
                if(!collection->value.IsNull() && !collection->value.Empty()) {
                    for(auto& element : collection->value.GetArray()) {
                        if(decode_file_url_by_key("url", url, IoDirection::IN, element)) {
                            if(input_feed_url.count(url) > 0) {
                                throw ConfigurationError("URL " + string(url) + " is used for both input and output");
                            }
                        }
                    }
                }
            }

        }
    }
};
void Environment::load_pheniqs_pg() {
    // ks_put_string(interface.name().c_str(), interface.name().size(), &_pheniqs_pg.ID);
    // ks_put_string(interface.name().c_str(), interface.name().size(), &_pheniqs_pg.PN);
    // ks_put_string(interface.application_version.c_str(), interface.application_version.size(), &_pheniqs_pg.VN);
    // ks_put_string(interface.full_command.c_str(), interface.full_command.size(), &_pheniqs_pg.CL);
};
void Environment::print_help(ostream& o) const {
    interface.print_help(o);
};
void Environment::print_version(ostream& o) const {
    interface.print_version(o);
};
void Environment::print_linted_instruction(ostream& o) const {
    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    _instruction.Accept(writer);
    o << buffer.GetString() << endl;
};
void Environment::print_instruction_validation(ostream& o) const {
    o << fixed << setprecision(19);
    o << "Environment " << endl << endl;
    o << "    Version                                     " << interface.application_version << endl;

    Decoder decoder;
    decode_value_by_key< Decoder >("decoder", decoder, _instruction);
    o << "    Multiplex barcode decoder                   " << decoder << endl;

    Platform platform;
    decode_value_by_key< Platform >("platform", platform, _instruction);
    o << "    Platform                                    " << platform << endl;

    bool long_read;
    decode_value_by_key< bool >("long read", long_read, _instruction);
    o << "    Long read                                   " << (long_read ? "enabled" : "disabled") << endl;

    bool disable_quality_control;
    decode_value_by_key< bool >("disable quality control", disable_quality_control, _instruction);
    o << "    Quality tracking                            " << (disable_quality_control ? "disabled" : "enabled") << endl;

    bool include_filtered;
    decode_value_by_key< bool >("include filtered", include_filtered, _instruction);
    o << "    Include non PF reads                        " << (include_filtered ? "enabled" : "disabled") << endl;

    uint8_t input_phred_offset;
    decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, _instruction);
    o << "    Input Phred offset                          " << to_string(input_phred_offset) << endl;

    uint8_t output_phred_offset;
    decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, _instruction);
    o << "    Output Phred offset                         " << to_string(output_phred_offset) << endl;

    uint64_t leading_segment_index;
    decode_value_by_key< uint64_t >("leading segment index", leading_segment_index, _instruction);
    o << "    Leading template segment                    " << to_string(leading_segment_index) << endl;

    // o << "    Undetermined reads                          " << (undetermined != NULL ? undetermined->alias() : "ignored") << endl;

    uint64_t multiplex_segment_cardinality;
    if(decode_value_by_key< uint64_t >("multiplex segment cardinality", multiplex_segment_cardinality, _instruction)) {
        if(multiplex_segment_cardinality) {
            uint64_t concatenated_multiplex_barcode_length;
            decode_value_by_key< uint64_t >("concatenated multiplex barcode length", concatenated_multiplex_barcode_length, _instruction);
            o << "    Multiplex barcoding cycles                  " << concatenated_multiplex_barcode_length << endl;

            vector< uint64_t > multiplex_barcode_length;
            decode_value_by_key< vector< uint64_t > >("multiplex barcode length", multiplex_barcode_length, _instruction);
            o << "    Multiplex barcode length                    ";
            for(const auto& value : multiplex_barcode_length) {
                o << value << ' ';
            } o << endl;

            o << "    Shortest multiplex distance                 " << _multiplex_barcode_distance.minimum_distance() << endl;

            if(decoder == Decoder::MDD || decoder == Decoder::BENCHMARK) {
                o << "    Tolerable multiplex shannon distance        ";
                for(const auto& value : _multiplex_barcode_set_distance) {
                    o << value.shannon_bound() << ' ';
                } o << endl;

                vector< uint8_t > multiplex_barcode_tolerance;
                decode_value_by_key< vector< uint8_t > >("multiplex barcode tolerance", multiplex_barcode_tolerance, _instruction);

                o << "    Multiplex barcode tolerance                 ";
                for(const auto& value : multiplex_barcode_tolerance) {
                    o << uint64_t(value) << ' ';
                } o << endl;

                uint8_t masking_threshold;
                if(decode_value_by_key< uint8_t >("masking_threshold", masking_threshold, _instruction)) {
                    o << "    Multiplex masking threshold                 " << (masking_threshold > 0 ? to_string(masking_threshold) : "NONE") << endl;
                }
            }

            if(decoder == Decoder::PAMLD || decoder == Decoder::BENCHMARK) {
                double multiplex_noise;
                decode_value_by_key< double >("multiplex noise", multiplex_noise, _instruction);
                o << "    Prior multiplex noise                       " << to_string(multiplex_noise) << endl;

                double multiplex_confidence;
                decode_value_by_key< double >("multiplex confidence", multiplex_confidence, _instruction);
                o << "    Multiplex confidence                        " << to_string(multiplex_confidence) << endl;

                double random_multiplex_barcode_probability;
                decode_value_by_key< double >("random multiplex barcode probability", random_multiplex_barcode_probability, _instruction);
                o << "    Random word frequency                       " << random_multiplex_barcode_probability << endl;

                double adjusted_multiplex_noise_probability;
                decode_value_by_key< double >("adjusted multiplex noise probability", adjusted_multiplex_noise_probability, _instruction);
                o << "    Adjusted multiplex noise probability        " << adjusted_multiplex_noise_probability << endl;
            }
        }
    }

    uint64_t molecular_segment_cardinality;
    if(decode_value_by_key< uint64_t >("molecular segment cardinality", molecular_segment_cardinality, _instruction)) {
        if(molecular_segment_cardinality) {
            uint64_t concatenated_molecular_barcode_length;
            decode_value_by_key< uint64_t >("concatenated molecular barcode length", concatenated_molecular_barcode_length, _instruction);
            o << "    Molecular barcoding cycles                  " << concatenated_molecular_barcode_length << endl;

            vector< uint64_t > molecular_barcode_length;
            decode_value_by_key< vector< uint64_t > >("molecular barcode length", molecular_barcode_length, _instruction);
            o << "    Molecular barcode length                    ";
            for(const auto& value : molecular_barcode_length) {
                o << value << ' ';
            } o << endl;
        }
    }

    int32_t buffer_capacity;
    decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, _instruction);
    o << "    Feed buffer capacity                        " << to_string(buffer_capacity) << endl;

    int32_t threads;
    decode_value_by_key< int32_t >("threads", threads, _instruction);
    o << "    Threads                                     " << to_string(threads) << endl;
    o << endl;

    o << "Transformation " << endl << endl;
    vector< Token > token_array;
    if(decode_value_by_key< vector< Token > >("token", token_array, _instruction)) {
        for(auto& token : token_array) {
            o << "    Token No." << token.index << endl;
            o << "        Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
            o << "        Pattern       " << string(token) << endl;
            o << "        Description   ";
            o << token.description() << endl;
            o << endl;
        }
    }

    list< Transform > template_transforms;
    if(decode_transform_array_by_key("template", template_transforms, _instruction, token_array)) {
        o << "    Template transform" << endl;
        for(const auto& transform : template_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }

    list< Transform > multiplex_barcode_transforms;
    if(decode_transform_array_by_key("template", multiplex_barcode_transforms, _instruction, token_array)) {
        o << "    Multiplex barcode transform" << endl;
        for(const auto& transform : multiplex_barcode_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }

    list< Transform > molecular_barcode_transforms;
    if(decode_transform_array_by_key("template", molecular_barcode_transforms, _instruction, token_array)) {
        o << "    Molecular barcode transform" << endl;
        for(const auto& transform : molecular_barcode_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }

    list< URL > input_url_array;
    if(decode_file_url_list_by_key("input", input_url_array, _instruction, IoDirection::IN)) {
        o << "Input " << endl << endl;
        uint64_t url_index(0);
        for(auto& url : input_url_array) {
            o << "    Input segment No." << url_index << " : " << url << endl;
            url_index++;
        }
        o << endl;
    }

    list< ChannelSpecification > channel_specifications;
    if(decode_value_by_key< list< ChannelSpecification > >("channel", channel_specifications, _instruction)) {
        o << "Channel " << endl << endl;
        for(const auto& specification : channel_specifications) {
            specification.describe(o);
        }
    }

    if(_display_distance && !_multiplex_barcode_distance.empty()) {
        o << "Barcode distance distribution" << endl << endl;
        for(const auto& metric : _multiplex_barcode_set_distance) {
            metric.describe(o);
        }

        if(_multiplex_barcode_set_distance.size() > 1) {
            _multiplex_barcode_distance.describe(o);
        }
    }
};
