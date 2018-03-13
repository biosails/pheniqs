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

#ifndef PHENIQS_ENVIRONMENT_H
#define PHENIQS_ENVIRONMENT_H

#include <set>
#include <unordered_map>

#include <htslib/kstring.h>
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/error/en.h>

#include "version.h"
#include "error.h"
#include "interface.h"
#include "model.h"
#include "feed.h"

using std::set;
using std::hash;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::pair;
using std::fixed;
using std::size_t;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::unordered_map;
using std::numeric_limits;
using std::istreambuf_iterator;

using rapidjson::Document;
using rapidjson::Value;
using rapidjson::SizeType;
using rapidjson::StringBuffer;
using rapidjson::PrettyWriter;


/*  Environment
*/
class Environment {
    public:
        CommandLine* interface;
        ProgramAction action;
        ProgramState state;
        string version;
        URL configuration_url;
        Decoder decoder;
        Platform platform;
        bool help_only;
        bool version_only;
        bool disable_quality_control;
        bool long_read;
        bool validate_only;
        bool lint_only;
        bool display_distance;
        bool include_filtered;
        string facility;
        string platform_model;
        string production_date;
        URL working_directory;
        URL base_input_url;
        URL base_output_url;
        string insert_size;
        uint32_t input_phred_offset;
        uint32_t output_phred_offset;
        uint32_t masking_threshold;
        uint32_t leading_segment_index;
        uint32_t threads;
        uint32_t transforms;
        uint32_t buffer_capacity;
        size_t total_multiplex_barcode_segments;
        size_t total_molecular_barcode_segments;
        size_t total_output_segments;
        size_t total_input_segments;
        double confidence;
        double noise;

        // input channel specification
        vector< URL > input_urls;
        URL input_url;

        // token specification
        vector< Token > tokens;

        // template transform specification
        vector< Transform > template_transforms;

        // multiplex and molecular barcode transform specification
        vector< Transform > multiplex_barcode_transforms;
        vector< Transform > molecular_barcode_transforms;
        vector< size_t > multiplex_barcode_length;
        vector< size_t > molecular_barcode_length;
        vector< uint8_t > multiplex_barcode_tolerance;

        // read group specification by id
        unordered_map< string, HeadRGAtom* > read_group_by_id;

        // feed specification by URL
        unordered_map< URL, FeedSpecification* > input_feed_specification_by_url;
        unordered_map< URL, FeedSpecification* > output_feed_specification_by_url;

        // channel specification
        vector< ChannelSpecification* > channel_specifications;

        // undetermined channel reference
        ChannelSpecification* undetermined;

        // computed
        double adjusted_noise_probability  = 0;
        double random_word_probability = 0;
        HeadPGAtom pg;

        Environment(int argc, char** argv);
        ~Environment();
        void print_help(ostream& o) const;
        void print_version(ostream& o) const;
        void print_configuration(ostream& o) const;
        void describe(ostream& o);
        inline size_t number_of_barcodes() const {
            return multiplex_barcode_distance.height();
        };
        inline size_t multiplex_barcode_width() const {
            return multiplex_barcode_distance.width();
        };
        inline size_t minimum_barcode_distance() const {
            return multiplex_barcode_distance.minimum_distance();
        };
        void validate_urls();
        void load_input_specification();
        void load_transformation();
        void load_channels();
        void calibrate(const URL& url);
        void encode(Document& document) const;
        ChannelSpecification* load_channel_from_rg(const HeadRGAtom& rg);
        FeedSpecification* discover_feed(const URL& url, const IoDirection& direction);

    private:
        vector< string > token_patterns;
        vector< string > template_patterns;
        vector< string > multiplex_barcode_patterns;
        vector< string > molecular_barcode_patterns;
        vector< Distance > multiplex_barcode_set_distance;
        Distance multiplex_barcode_distance;

        void load(int argc, char** argv);
        void load_configuration_file(const URL& url);
        void load_url_node(const Value& node, const IoDirection& direction, URL& url);
        void load_token_node(const Value& node);
        void load_transform_node(const Value& node, vector< string >& container);
        void load_barcode_node(const Value& node, Barcode& barcode);
        void load_channel_node(const Value& node);
        void load_read_group_node(const Value& node);
        void load_defaults();
        void load_urls();
        void load_token(const string& pattern);
        void load_transform(const string& pattern, vector< Transform >& container, const size_t& index);
        void load_barcode_tolerance();
        void load_thread_model();
        void load_undetermined();
        void load_multiplex_barcodes();
        void load_prior();
        void load_output_specification();
        void validate_io();
        void probe();
        void validate();

        inline void set_state(const ProgramState code) {
            state = code;
        };
        inline void set_selected_action(string& name) {
            if (!name.compare("demux")) {
                action = ProgramAction::DEMULTIPLEX;
            } else if (!name.compare("quality")) {
                action = ProgramAction::QUALITY;
            }
        };
        inline void set_help_only(const bool* value) {
            if(value!=NULL) {
                help_only = *value;
            }
        };
        inline void set_version_only(const bool* value) {
            if(value!=NULL) {
                version_only = *value;
            }
        };
        inline void set_long_read(const bool* value) {
            if(value!=NULL) {
                long_read = *value;
            }
        };
        inline void set_disable_quality_control(const bool* value) {
            if(value!=NULL) {
                disable_quality_control = *value;
            }
        };
        inline void set_validate_only(const bool* value) {
            if(value!=NULL) {
                validate_only = *value;
            }
        };
        inline void set_lint_only(const bool* value) {
            if(value!=NULL) {
                lint_only = *value;
            }
        };
        inline void set_display_distance(const bool* value) {
            if(value!=NULL) {
                display_distance = *value;
            }
        };
        inline void set_include_filtered(const bool* value) {
            if(value!=NULL) {
                include_filtered = *value;
            }
        };
        inline void set_decoder(const string* value) {
            if(value!=NULL) {
                set_decoder(value->c_str());
            }
        };
        inline void set_decoder(const char* value) {
            if(value!=NULL) {
                value >> decoder;
                if(decoder == Decoder::UNKNOWN) {
                    throw ConfigurationError("bad decoder value " + string(value));
                }
            }
        };
        inline void set_platform(const string* value) {
            if(value!=NULL) {
                set_platform(value->c_str());
            }
        };
        inline void set_platform(const char* value) {
            if(value!=NULL) {
                value >> platform;
                if(platform == Platform::UNKNOWN) {
                    throw ConfigurationError("bad platform value " + string(value));
                }
            }
        };
        inline void set_configuration_path(const string* value) {
            if(value!=NULL) {
                configuration_url.parse(*value, IoDirection::IN);
            }
        };
        inline void set_base_input_path(const string* value) {
            if(value!=NULL) {
                base_input_url.set_directory(*value);
            }
        };
        inline void set_base_output_path(const string* value) {
            if(value!=NULL) {
                base_output_url.set_directory(*value);
            }
        };
        inline void set_threads(const long* value) {
            if(value!=NULL) {
                threads = size_t(*value);
            }
        };
        inline void set_pivot_threads(const long* value) {
            if(value!=NULL) {
                transforms = size_t(*value);
            }
        };
        inline void set_buffer_capacity(const long* value) {
            if(value!=NULL) {
                buffer_capacity = size_t(*value);
            }
        };
        inline void set_noise(const double* value) {
            if(value!=NULL) {
                noise = *value;
            }
        };
        inline void set_confidence(const double* value) {
            if(value!=NULL) {
                confidence = *value;
            }
        };
        inline void set_leading_segment_index(const long* value) {
            if(value!=NULL) {
                leading_segment_index = uint8_t(*value);
            }
        };
        inline void set_decoder_masking_threshold(const long* value) {
            if(value!=NULL) {
                masking_threshold = uint8_t(*value);
            }
        };
        inline void set_input_path(const string* value) {
            if(value!=NULL) {
                input_url.parse(*value, IoDirection::IN);
            }
        };
        inline void decode_directory_node(const Value& node, const Value::Ch* name, URL& value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsString()) {
                    value.set_directory(string(element->value.GetString(), element->value.GetStringLength()));
                } else { throw ConfigurationError(string(name) + " element must be a string"); }
            }
        };
        inline void decode_string_node(const Value& node, const Value::Ch* name, string& value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsString()) {
                    value.assign(element->value.GetString(), element->value.GetStringLength());
                } else { throw ConfigurationError(string(name) + " element must be a string"); }
            }
        };
        inline void decode_string_node(const Value& node, const Value::Ch* name, kstring_t* value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsString()) {
                    kputsn(element->value.GetString(), element->value.GetStringLength(), value);
                } else { throw ConfigurationError(string(name) + " element must be a string"); }
            }
        };
        inline void decode_uint_node(const Value& node, const Value::Ch* name, uint32_t& value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsUint()) {
                    value = element->value.GetUint();
                } else { throw ConfigurationError(string(name) + " element must be an unsigned integer"); }
            }
        };
        inline void decode_bool_node(const Value& node, const Value::Ch* name, bool& value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsBool()) {
                    value = element->value.GetBool();
                } else { throw ConfigurationError(string(name) + " element must be a boolean"); }
            }
        };
        inline void decode_double_node(const Value& node, const Value::Ch* name, double& value) {
            Value::ConstMemberIterator element = node.FindMember(name);
            if (element != node.MemberEnd()) {
                if(element->value.IsNumber()) {
                    value = element->value.GetDouble();
                } else { throw ConfigurationError(string(name) + " element must be numeric"); }
            }
        };
        inline void decode_read_group(HeadRGAtom& rg, const Value& node, const Value::Ch* key) {
            if (node.IsObject()) {
                decode_string_node(node, key, &rg.ID);
                decode_string_node(node, "PI", &rg.PI);
                decode_string_node(node, "LB", &rg.LB);
                decode_string_node(node, "SM", &rg.SM);
                decode_string_node(node, "PU", &rg.PU);
                decode_string_node(node, "CN", &rg.CN);
                decode_string_node(node, "DS", &rg.DS);
                decode_string_node(node, "DT", &rg.DT);
                decode_string_node(node, "PL", &rg.PL);
                decode_string_node(node, "PM", &rg.PM);
                decode_string_node(node, "PG", &rg.PG);
                decode_string_node(node, "FO", &rg.FO);
                decode_string_node(node, "KS", &rg.KS);
            }
        };
        inline void encode_transform(Document& document, Value& node, const vector< Transform >& container, const string& key) const {
            if(!container.empty()) {
                Document::AllocatorType& allocator = document.GetAllocator();

                Value collection;
                collection.SetArray();

                size_t index = 0;
                string current;
                for(auto& transform : container) {
                    if(transform.output_segment_index != index) {
                        collection.PushBack(Value(current.c_str(), current.size(), allocator).Move(), allocator);
                        current.clear();
                        index++;
                    }
                    if(current.size() > 0) {
                        current.push_back(':');
                    }
                    current.append(string(transform));
                }
                collection.PushBack(Value(current.c_str(), current.size(), allocator).Move(), allocator);
                node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), collection, allocator);
            }
        };
};
#endif /* PHENIQS_ENVIRONMENT_H */