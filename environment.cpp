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
#include "configuration.h"

static inline string getcwd() {
    char* buffer;
    char* temp;
    string cwd;
    size_t size = 128;

    if((buffer = (char*)malloc(size)) == NULL) {
        throw InternalError("out of memory");
    }
    while(getcwd(buffer, size) == NULL) {
        switch(errno) {
            case EACCES:
                free(buffer);
                throw IOError("insufficient permission to probe working directory");
                break;

            case ERANGE:
                size *= 2;
                if((temp = (char*)realloc(buffer, size)) == NULL) {
                    free(buffer);
                    throw InternalError("out of memory");
                } else {
                    buffer = temp;
                }
                break;

            default:
                throw InternalError("error " + to_string(errno) + " when probing working directory");
                break;
        }
    }
    cwd.assign(buffer);
    free(buffer);
    return cwd;
}

/*  Environment
*/
Environment::Environment(int argc, char** argv) :
    state(ProgramState::OK),
    decoder(Decoder::UNKNOWN),
    platform(Platform::UNKNOWN),
    help_only(false),
    version_only(false),
    disable_quality_control(false),
    long_read(false),
    validate_only(false),
    lint_only(false),
    display_distance(false),
    include_filtered(false),
    input_phred_offset(numeric_limits<uint32_t>::max()),
    output_phred_offset(numeric_limits<uint32_t>::max()),
    masking_threshold(0),
    leading_segment_index(0),
    threads(1),
    transforms(numeric_limits<uint32_t>::max()),
    buffer_capacity(numeric_limits<uint32_t>::max()),
    total_multiplex_barcode_segments(0),
    total_molecular_barcode_segments(0),
    total_output_segments(0),
    total_input_segments(0),
    confidence(0),
    noise(0),
    undetermined(NULL) {

    load(argc, argv);
};
Environment::~Environment() {
    for(const auto& record : input_feed_specification_by_url) {
        delete record.second;
    }
    input_feed_specification_by_url.clear();

    for(const auto& record : output_feed_specification_by_url) {
        delete record.second;
    }
    output_feed_specification_by_url.clear();

    for(auto specification : channel_specifications) {
        delete specification;
    }
    channel_specifications.clear();

    delete interface;
};
void Environment::load(int argc, char** argv) {
    try {
        version = PHENIQS_VERSION;
        const string content(reinterpret_cast<const char *>(interface_json), interface_json_len);

        try {
            interface = new CommandLine(content.c_str(), version);
        } catch(ConfigurationError& e) {
            set_state(ProgramState::INTERNAL_ERROR);
            throw e;
        }

        string cwd = getcwd();
        working_directory.set_directory(string(cwd.c_str(), cwd.size()));

        interface->load(argc, argv);

        if(interface->help_triggered()) {
            set_state(ProgramState::HELP);

        } else if(interface->version_triggered()) {
            set_state(ProgramState::VERSION);

        } else {

            // set the selected action
            set_selected_action(interface->get_selected_action());

            // first get the configuration file URL and load it
            set_configuration_path(interface->get_string("configuration path"));
            load_configuration_file(configuration_url);

            // then decode the rest of the parameters to override whatever came from the configuration file
            set_help_only(interface->get_boolean("help"));
            set_version_only(interface->get_boolean("version"));
            set_disable_quality_control(interface->get_boolean("quality"));
            set_long_read(interface->get_boolean("long read"));
            set_validate_only(interface->get_boolean("validate"));
            set_lint_only(interface->get_boolean("lint"));
            set_display_distance(interface->get_boolean("display distance"));
            set_include_filtered(interface->get_boolean("filtered"));
            set_decoder(interface->get_string("decoder"));
            set_platform(interface->get_string("platform"));
            set_base_input_path(interface->get_string("base input path"));
            set_base_output_path(interface->get_string("base output path"));
            set_threads(interface->get_integer("threads"));
            set_pivot_threads(interface->get_integer("transforms"));
            set_buffer_capacity(interface->get_integer("buffer capacity"));
            set_confidence(interface->get_decimal("confidence"));
            set_noise(interface->get_decimal("noise"));
            set_leading_segment_index(interface->get_integer("leading segment index"));
            set_decoder_masking_threshold(interface->get_integer("masking threshold"));
            set_input_path(interface->get_string("input path"));

            if(state == ProgramState::OK) {
                load_defaults();
                load_urls();
                load_thread_model();

                switch (action) {
                    case ProgramAction::DEMULTIPLEX: {
                        load_transformation();
                        load_barcode_tolerance();
                        load_undetermined();
                        load_multiplex_barcodes();
                        load_prior();
                        validate();
                        load_channels();
                        load_input_specification();
                        load_output_specification();
                        validate_io();
                        probe();
                        break;
                    };

                    case ProgramAction::QUALITY: {
                        break;
                    };

                    default: break;
                }
            }
        }

    } catch(IOError& e) {
        set_state(ProgramState::IO_ERROR);
        throw e;

    } catch(CommandLineError& e) {
        set_state(ProgramState::COMMAND_LINE_ERROR);
        throw e;

    } catch(ConfigurationError& e) {
        set_state(ProgramState::CONFIGURATION_ERROR);
        throw e;
    }
};
void Environment::print_help(ostream& o) const {
    interface->print_help(o);
};
void Environment::print_version(ostream& o) const {
    interface->print_version(o);
};
void Environment::print_configuration(ostream& o) const {
    Document document;
    encode(document);
    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    document.Accept(writer);
    o << buffer.GetString() << endl;
};
void Environment::describe(ostream& o) {
    o << fixed << setprecision(19);
    o << "Environment " << endl << endl;
    {
        o << "    Version                                     " << version << endl;
        o << "    Barcode decoder                             " << decoder << endl;
        o << "    Platform                                    " << platform << endl;
        o << "    Long read                                   " << (long_read ? "enabled" : "disabled") << endl;
        o << "    Quality tracking                            " << (disable_quality_control ? "disabled" : "enabled") << endl;
        o << "    Leading template segment                    " << to_string(leading_segment_index) << endl;
        o << "    Input Phred offset                          " << to_string(input_phred_offset) << endl;
        o << "    Output Phred offset                         " << to_string(output_phred_offset) << endl;
        o << "    Include non PF reads                        " << (include_filtered ? "enabled" : "disabled") << endl;
        o << "    Undetermined reads                          " << (undetermined != NULL ? undetermined->alias() : "ignored") << endl;
    }
    if (total_multiplex_barcode_segments > 0) {
        o << "    Multiplex barcoding cycles                  " << multiplex_barcode_width() << endl;
        o << "    Shortest multiplex distance                 " << minimum_barcode_distance() << endl;
        o << "    Multiplex barcode length                    ";
            for(const auto& length : multiplex_barcode_length) {
                o << size_t(length) << ' ';
            } o << endl;

        if(decoder == Decoder::MDD || decoder == Decoder::BENCHMARK) {
        o << "    Tolerable shannon distance                  ";
            for(const auto metric : multiplex_barcode_set_distance) {
                o << metric.shannon_bound() << ' ';
            } o << endl;
        o << "    Barcode distance tolerance                  ";
            for(const auto& tolerance : multiplex_barcode_tolerance) {
                o << size_t(tolerance) << ' ';
            } o << endl;
        o << "    Multiplex masking threshold                 " << (masking_threshold > 0 ? to_string(masking_threshold) : "NONE") << endl;
        }

        if(decoder == Decoder::PAMLD || decoder == Decoder::BENCHMARK) {
        o << "    Prior noise frequency                       " << to_string(noise) << endl;
        o << "    Multiplex confidence                        " << to_string(confidence) << endl;
        o << "    Random word frequency                       " << random_word_probability << endl;
        }
    }
    if (total_molecular_barcode_segments > 0) {
        o << "    Molecular barcode length                    ";
            for(const auto& length : molecular_barcode_length) {
                o << size_t(length) << ' ';
            } o << endl;
    }
    {
        o << "    Feed buffer capacity                        " << to_string(buffer_capacity) << endl;
        o << "    Transforming threads                        " << to_string(transforms) << endl;
        o << "    Threads                                     " << to_string(threads) << endl;
    }
    o << endl;

    o << "Transformation " << endl << endl;
    if (!tokens.empty()) {
        for(const auto& token: tokens) {
            o << "    Token No." << token.index << endl;
            o << "        Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
            o << "        Pattern       " << string(token) << endl;
            o << "        Description   ";
            o << token.description() << endl;
            o << endl;
        }
    }
    if (!template_transforms.empty()) {
        o << "    Template transform" << endl;
        for(const auto& transform : template_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }
    if (!multiplex_barcode_transforms.empty()) {
        o << "    Multiplex barcode transform" << endl;
        for(const auto& transform : multiplex_barcode_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }
    if (!molecular_barcode_transforms.empty()) {
        o << "    Molecular barcode transform" << endl;
        for(const auto& transform : molecular_barcode_transforms) {
            o << "        " << transform.description() << endl;
        }
        o << endl;
    }

    o << "Input " << endl << endl;
    if (!input_urls.empty()) {
        for (size_t i = 0; i < input_urls.size(); i++) {
            o << "    Input segment No." << i << " : " << input_urls[i] << endl;
        }
        o << endl;
    }

    o << "Channel " << endl << endl;
    if(!channel_specifications.empty()) {
        for(const auto specification : channel_specifications) {
            specification->describe(o);
        }
    }

    o << "Feed " << endl << endl;
    for(const auto& item : input_feed_specification_by_url) {
        item.second->describe(o);
    }
    for(const auto& item : output_feed_specification_by_url) {
        item.second->describe(o);
    }
    if(display_distance && !multiplex_barcode_distance.empty()) {
        o << "Barcode distance distribution" << endl << endl;
        for(const auto metric : multiplex_barcode_set_distance) {
            metric.describe(o);
        }

        if(multiplex_barcode_set_distance.size() > 1) {
            multiplex_barcode_distance.describe(o);
        }
    }
};
void Environment::encode(Document& document) const {
    document.SetObject();
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;
    Value collection;

    if(!facility.empty()) {
        v.SetString(facility.c_str(), facility.size(), allocator);
        document.AddMember("CN", v, allocator);
    }
    if(!platform_model.empty()) {
        v.SetString(platform_model.c_str(), platform_model.size(), allocator);
        document.AddMember("PM", v, allocator);
    }
    if(!production_date.empty()) {
        v.SetString(production_date.c_str(), production_date.size(), allocator);
        document.AddMember("DT", v, allocator);
    }
    if(!insert_size.empty()) {
        v.SetString(insert_size.c_str(), insert_size.size(), allocator);
        document.AddMember("PI", v, allocator);
    }
    if(!base_input_url.empty()) {
        // base_input_url.encode(document, v);
        document.AddMember("base input path", v.Move(), allocator);
    }
    if(!base_output_url.empty()) {
        // base_output_url.encode(document, v);
        document.AddMember("base output path", v.Move(), allocator);
    }
    if(platform != Platform::UNKNOWN) {
        string buffer;
        buffer << platform;
        v.SetString(buffer.c_str(), buffer.size(), allocator);
        document.AddMember("PL", v, allocator);
    }
    if(decoder != Decoder::UNKNOWN) {
        string buffer;
        buffer << decoder;
        v.SetString(buffer.c_str(), allocator);
        document.AddMember("decoder", v, allocator);
    }

    v.SetBool(disable_quality_control);
    document.AddMember("disable quality control", v, allocator);

    v.SetBool(long_read);
    document.AddMember("long read", v, allocator);

    v.SetBool(include_filtered);
    document.AddMember("include filtered", v, allocator);

    v.SetUint64(input_phred_offset);
    document.AddMember("input phred offset", v, allocator);

    v.SetUint64(output_phred_offset);
    document.AddMember("output phred offset", v, allocator);

    v.SetUint64(masking_threshold);
    document.AddMember("decoder masking threshold", v, allocator);

    v.SetUint64(leading_segment_index);
    document.AddMember("leading segment", v, allocator);

    v.SetUint64(threads);
    document.AddMember("threads", v, allocator);

    v.SetUint64(transforms);
    document.AddMember("transforms", v, allocator);

    v.SetUint64(buffer_capacity);
    document.AddMember("buffer capacity", v, allocator);

    v.SetDouble(confidence);
    document.AddMember("confidence", v, allocator);

    v.SetDouble(noise);
    document.AddMember("noise", v, allocator);

    if(!input_urls.empty()) {
        collection.SetArray();
        for(auto& url : input_urls) {
            // url.encode(document, v);
            collection.PushBack(v, allocator);
        }
        document.AddMember("input", collection, allocator);
    }
    if(!multiplex_barcode_tolerance.empty()) {
        collection.SetArray();
        for(auto& tolerance : multiplex_barcode_tolerance) {
            v.SetUint64(tolerance);
            collection.PushBack(v, allocator);
        }
        document.AddMember("distance tolerance", collection, allocator);
    }
    if(!tokens.empty()) {
        collection.SetArray();
        for(auto& token : tokens) {
            string buffer(token);
            v.SetString(buffer.c_str(), buffer.size(), allocator);
            collection.PushBack(v, allocator);
        }
        document.AddMember("token", collection, allocator);
    }
    encode_transform(document, document, template_transforms, "template");
    encode_transform(document, document, multiplex_barcode_transforms, "multiplex barcode");
    encode_transform(document, document, molecular_barcode_transforms, "molecular barcode");

    if(!channel_specifications.empty()) {
        collection.SetArray();
        for(auto specification : channel_specifications) {
            specification->encode(document, collection);
        }
        document.AddMember("channel", collection, allocator);
    }
};
void Environment::validate_urls() {
    for (auto& url : input_urls) {
        if(!url.is_readable()) {
            set_state(ProgramState::IO_ERROR);
            throw IOError("could not open " + string(url) + " for reading");
        }
    }

    for(const auto specification : channel_specifications) {
        for (const auto& url : specification->output_urls) {
            if(!url.is_writable()) {
                set_state(ProgramState::IO_ERROR);
                throw IOError("could not open " + string(url) + " for writing");
            }
        }
    }
};

/*  configuration file parsing
*/
void Environment::load_configuration_file(const URL& url) {
    Document document;
    if (!url.empty()) {
        if(url.is_readable()) {
            ifstream file(url.path());
            string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
            file.close();
            if (!document.Parse(content.c_str()).HasParseError()) {
                if (document.IsObject()) {
                    Value::ConstMemberIterator element;
                    decode_string_node(document, "CN", facility);
                    decode_string_node(document, "PM", platform_model);
                    decode_string_node(document, "DT", production_date);
                    decode_string_node(document, "PI", insert_size);
                    decode_directory_node(document, "base input path", base_input_url);
                    decode_directory_node(document, "base output path", base_output_url);
                    decode_uint_node(document, "input phred offset", input_phred_offset);
                    decode_uint_node(document, "output phred offset", output_phred_offset);
                    decode_uint_node(document, "decoder masking threshold", masking_threshold);
                    decode_uint_node(document, "leading segment", leading_segment_index);
                    decode_bool_node(document, "disable quality control", disable_quality_control);
                    decode_bool_node(document, "long read", long_read);
                    decode_bool_node(document, "include filtered", include_filtered);
                    decode_bool_node(document, "validate only", validate_only);
                    decode_bool_node(document, "lint only", lint_only);
                    decode_bool_node(document, "display distance", display_distance);
                    decode_uint_node(document, "threads", threads);
                    decode_uint_node(document, "transforms", transforms);
                    decode_uint_node(document, "buffer capacity", buffer_capacity);
                    decode_double_node(document, "confidence", confidence);
                    decode_double_node(document, "noise", noise);

                    element = document.FindMember("decoder");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsString()) {
                            set_decoder(element->value.GetString());
                        } else { throw ConfigurationError("decoder element must be a string"); }
                    }
                    element = document.FindMember("PL");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsString()) {
                            set_platform(element->value.GetString());
                        } else { throw ConfigurationError("platform element must be a string"); }
                    }
                    element = document.FindMember("input");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsArray()) {
                            for (SizeType i = 0; i < element->value.Size(); i++) {
                                URL url;
                                load_url_node(element->value[i], IoDirection::IN, url);
                                input_urls.push_back(url);
                            }
                        } else { throw ConfigurationError("input element must be an array"); }
                    }

                    element = document.FindMember("token");
                    if (element != document.MemberEnd()) {
                        load_token_node(element->value);
                    }
                    element = document.FindMember("template");
                    if (element != document.MemberEnd()) {
                        load_transform_node(element->value, template_patterns);
                    }
                    element = document.FindMember("multiplex barcode");
                    if (element != document.MemberEnd()) {
                        load_transform_node(element->value, multiplex_barcode_patterns);
                    }
                    element = document.FindMember("molecular barcode");
                    if (element != document.MemberEnd()) {
                        load_transform_node(element->value, molecular_barcode_patterns);
                    }
                    element = document.FindMember("read group");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsArray()) {
                            for (SizeType i = 0; i < element->value.Size(); i++) {
                                try {
                                    load_read_group_node(element->value[i]);
                                } catch(ConfigurationError& e) {
                                    e.message = "read group in position " + to_string(i) + " " + e.message;
                                    throw e;
                                }
                            }
                        } else { throw ConfigurationError("read group element must be an array"); }
                    }
                    element = document.FindMember("channel");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsArray()) {
                            for (SizeType i = 0; i < element->value.Size(); i++) {
                                try {
                                    load_channel_node(element->value[i]);
                                } catch(ConfigurationError& e) {
                                    e.message = "channel in position " + to_string(i) + " " + e.message;
                                    throw e;
                                }
                            }
                        } else { throw ConfigurationError("channel element must be an array"); }
                    }
                    element = document.FindMember("distance tolerance");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsArray()) {
                            for (SizeType i = 0; i < element->value.Size(); i++) {
                                if (element->value[i].IsUint()) {
                                    multiplex_barcode_tolerance.push_back(element->value[i].GetUint());
                                } else {
                                    throw ConfigurationError("distance tolerance at position " + to_string(i) + " must be a positive integer");
                                }
                            }
                        } else { throw ConfigurationError("distance tolerance element must be an array"); }
                    }
                } else {
                    throw ConfigurationError("configuration root node must be a dictionary");
                }
            } else {
                string message(GetParseError_En(document.GetParseError()));
                message += " at position ";
                message += to_string(document.GetErrorOffset());
                throw ConfigurationError(message);
            }
        } else {
            throw IOError("could not read configuration from " + string(url));
        }
    }
};
void Environment::load_url_node(const Value& node, const IoDirection& direction, URL& url) {
    url.clear();
    if(node.IsString() || node.IsObject()) {
        try {
            if (node.IsString()) {
                url.parse(string(node.GetString(), node.GetStringLength()), direction);

            } else {
                string buffer;
                decode_string_node(node, "path", buffer);
                if(!buffer.empty()) {
                    url.parse(buffer, direction);
                } else {
                    throw ConfigurationError("URL element must contain a non empty path element");
                }

                buffer.clear();
                decode_string_node(node, "type", buffer);
                if(!buffer.empty()) {
                    url.set_type(buffer.c_str());
                }

                buffer.clear();
                decode_string_node(node, "compression", buffer);
                if(!buffer.empty()) {
                    url.set_compression(buffer);
                }
            }
        } catch(ConfigurationError& e) {
            url.clear();
            throw e;
        }
    } else {
        throw ConfigurationError("URL element must be either a string or a dictionary");
    }
};
void Environment::load_token_node(const Value& node) {
    if (node.IsArray()) {
        for (SizeType i = 0; i < node.Size(); i++) {
            if (node[i].IsString()) {
                token_patterns.emplace_back(node[i].GetString(), node[i].GetStringLength());
            } else {
                throw ConfigurationError("token array element at position " + to_string(i) + " must be a string");
            }
        }
    } else { throw ConfigurationError("token element must be an array"); }
};
void Environment::load_transform_node(const Value& node, vector< string >& container) {
    if (node.IsArray()) {
        for (SizeType i = 0; i < node.Size(); i++) {
            if (node[i].IsString()) {
                container.emplace_back(node[i].GetString(), node[i].GetStringLength());
            } else {
                throw ConfigurationError("transform at position " + to_string(i) + " must be a string");
            }
        }
    } else { throw ConfigurationError("transform element must be an array"); }
};
void Environment::load_read_group_node(const Value& node) {
    if (node.IsObject()) {
        HeadRGAtom* rg = new HeadRGAtom();
        decode_read_group(*rg, node, "ID");
        if(rg->ID.l > 0) {
            string id(*rg);
            auto record = read_group_by_id.find(id);
            if (record == read_group_by_id.end()) {
                read_group_by_id.emplace(make_pair(id, rg));
            } else {
                record->second->expand(*rg);
                delete rg;
            }
        } else {
            delete rg;
            throw ConfigurationError("is missing an ID");
        }
    } else {
        throw ConfigurationError("must be a dictionary");
    }
};
void Environment::load_channel_node(const Value& node) {
    if (node.IsObject()) {
        ChannelSpecification* specification = new ChannelSpecification(channel_specifications.size());
        decode_read_group(specification->rg, node, "RG");
        decode_bool_node(node, "undetermined", specification->undetermined);
        decode_string_node(node, "FS", &specification->FS);
        decode_string_node(node, "CO", &specification->CO);
        decode_double_node(node, "concentration", specification->concentration);

        Value::ConstMemberIterator element;
        element = node.FindMember("barcode");
        if (element != node.MemberEnd()) {
            try {
                load_barcode_node(element->value, specification->multiplex_barcode);
            } catch(ConfigurationError& e) {
                delete specification;
                throw e;
            }
        }
        element = node.FindMember("output");
        if (element != node.MemberEnd()) {
            if(element->value.IsArray()) {
                for (SizeType i = 0; i < element->value.Size(); i++) {
                    URL url;
                    try {
                        load_url_node(element->value[i], IoDirection::OUT, url);
                    } catch(ConfigurationError& e) {
                        delete specification;
                        throw e;
                    }
                    specification->output_urls.emplace_back(url);
                }
            } else { 
                delete specification;
                throw ConfigurationError("output element must be an array");
            }
        }
        channel_specifications.push_back(specification);
    } else {
        throw ConfigurationError("must be a dictionary");
    }
};
void Environment::load_barcode_node(const Value& node, Barcode& barcode) {
    if(node.IsArray()) {
        for (SizeType i = 0; i < node.Size(); i++) {
            if (node[i].IsString()) {
                string buffer(node[i].GetString(), node[i].GetStringLength());
                for(auto& c : buffer) {
                    if(!is_iupac_strict(c)) {
                        throw ConfigurationError("contains an ambiguous nucleotide in barcode " + buffer);
                    }
                }
                barcode.fill(i, buffer.c_str(), buffer.size());
            } else {
                throw ConfigurationError("barcode segment " + to_string(i) + " must be a string");
            }
        }
    } else {
        throw ConfigurationError("barcode element must be an array");
    }
};
ChannelSpecification* Environment::load_channel_from_rg(const HeadRGAtom& rg) {
    ChannelSpecification* specification = new ChannelSpecification(channel_specifications.size());
    specification->rg = rg;
    channel_specifications.push_back(specification);
    return specification;
};
/*  loading data structures
*/
void Environment::load_defaults() {
    kputsn(interface->name().c_str(), interface->name().size(), &pg.ID);
    kputsn(interface->name().c_str(), interface->name().size(), &pg.PN);
    kputsn(interface->get_full_command().c_str(), interface->get_full_command().size(), &pg.CL);
    kputsn(version.c_str(), version.size(), &pg.VN);

    // Default decoder
    if(decoder == Decoder::UNKNOWN) {
        decoder = Decoder::PAMLD;
    }
    // Default platform
    if(platform == Platform::UNKNOWN) {
        platform = Platform::ILLUMINA;
    }

    // Default IO feeds cache capacity
    if(buffer_capacity == numeric_limits<uint32_t>::max()) {
        buffer_capacity = DEFAULT_BUFFER_CAPACITY;
    }

    // Default Phred decoding offset
    if(input_phred_offset == numeric_limits<uint32_t>::max()) {
        input_phred_offset = DEFAULT_PHRED_OFFSET;
    }

    if(output_phred_offset == numeric_limits<uint32_t>::max()) {
        output_phred_offset = DEFAULT_PHRED_OFFSET;
    }
};
void Environment::load_urls() {
    if(base_input_url.empty()) {
        base_input_url = working_directory;
    }
    for(auto& url : input_urls) {
        url.relocate(base_input_url);
    }

    if(base_output_url.empty()) {
        base_output_url = working_directory;
    }
    for(auto specification : channel_specifications) {
        for(auto& url : specification->output_urls) {
            url.relocate(base_output_url);
        }
    }
};
void Environment::load_transformation() {
    total_input_segments = input_urls.size();

    tokens.reserve(token_patterns.size());
    for(const auto& pattern : token_patterns) {
        load_token(pattern);
    }

    for(size_t i = 0; i < template_patterns.size(); i++) {
        load_transform(template_patterns[i], template_transforms, i);
    }
    total_output_segments = template_patterns.size();


    // multiplex barcode
    total_multiplex_barcode_segments = multiplex_barcode_patterns.size();
    for(size_t i = 0; i < total_multiplex_barcode_segments; i++) {
        load_transform(multiplex_barcode_patterns[i], multiplex_barcode_transforms, i);
    }
    multiplex_barcode_length.resize(total_multiplex_barcode_segments);
    for(auto& transform : multiplex_barcode_transforms) {
        if(transform.token.constant()) {
            if(!transform.token.empty()) {
                multiplex_barcode_length[transform.output_segment_index] += transform.token.length();
            } else {
                throw ConfigurationError("multiplex barcode token " + string(transform.token) + " is empty");
            }
        } else {
            throw ConfigurationError("multiplex barcode token " + string(transform.token) + " is not fixed width");
        }
    }

    // molecular barcode
    total_molecular_barcode_segments = molecular_barcode_patterns.size();
    for(size_t i = 0; i < total_molecular_barcode_segments; i++) {
        load_transform(molecular_barcode_patterns[i], molecular_barcode_transforms, i);
    }
    molecular_barcode_length.resize(total_molecular_barcode_segments);
    for(auto& transform : molecular_barcode_transforms) {
        if(transform.token.constant()) {
            if(!transform.token.empty()) {
                molecular_barcode_length[transform.output_segment_index] += transform.token.length();
            } else {
                throw ConfigurationError("molecular barcode token " + string(transform.token) + " is empty");
            }
        } else {
            throw ConfigurationError("molecular barcode token " + string(transform.token) + " is not fixed width");
        }
    }

    for(auto specification : channel_specifications) {
        if(!specification->output_urls.empty() && specification->output_urls.size() != total_output_segments) {
            if(specification->output_urls.size() == 1) {
                while(specification->output_urls.size() < total_output_segments) {
                    specification->output_urls.push_back(specification->output_urls[0]);
                }
            } else {
                throw ConfigurationError("incorrect number of output paths in channel " + to_string(specification->index));
            }
        }
    }
};
void Environment::load_token(const string& pattern) {
    size_t input_segment_index = numeric_limits<size_t>::max();
    int32_t start = numeric_limits<int32_t>::max();
    int32_t end = numeric_limits<int32_t>::max();
    bool end_terminated = true;

    bool sign = true;
    uint8_t literal = 0;
    size_t position = 0;
    int32_t value = numeric_limits<int32_t>::max();
    while(true) {
        const char& c = pattern[position];
        switch(c) {
            case ':':
            case '\0':
                switch(literal) {
                    case 0: {
                        if(value == numeric_limits<int32_t>::max()) {
                            throw ConfigurationError("token must explicitly specify an input segment reference");
                        } else if(value < 0) {
                            throw ConfigurationError("input segment reference must be a positive number");
                        } else {
                            input_segment_index = size_t(value);
                        }
                        break;
                    };
                    case 1:{
                        if(value == numeric_limits<int32_t>::max()) {
                            start = 0;
                        } else {
                            start = sign ? value : -value;
                        }
                        break;
                    };
                    case 2: {
                        if(value == numeric_limits<int32_t>::max()) {
                            end = 0;
                            end_terminated = false;
                        } else {
                            end = sign ? value : -value;
                        }
                        break;
                    };
                    default:
                        throw ConfigurationError("illegal token syntax " + pattern);
                        break;
                }
                value = numeric_limits<int32_t>::max();
                sign = true;
                literal++;
                break;
            case '-':
                sign = false;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9': {
                if(value == numeric_limits<int32_t>::max()) {
                    value = c - '0';
                } else {
                    value = value * 10 + (c - '0');
                }
                break;
            };
            default:
                throw ConfigurationError("illegal character " + to_string(c) + " in token");
                break;
        }
        if(c == '\0') { break; }
        position++;
    }
    if(input_segment_index < total_input_segments) {
        tokens.emplace_back(tokens.size(), input_segment_index, start, end, end_terminated);
    } else {
        throw ConfigurationError("invalid input segment reference " + to_string(input_segment_index) + " in token " + pattern);
    } 
};
void Environment::load_transform(const string& pattern, vector< Transform >& container, const size_t& index) {
    size_t position = 0;
    size_t v = numeric_limits<size_t>::max();
    LeftTokenOperator left = LeftTokenOperator::NONE;
    while(true) {
        const char& c = pattern[position];
        switch(c) {
            case ':':
            case '\0':
                if(v == numeric_limits<size_t>::max()) {
                    throw ConfigurationError("transform must explicitly specify a token reference");
                } else if(!(v < tokens.size())) {
                    throw ConfigurationError("invalid token reference " + to_string(v) + " in transform");
                } else {
                    const Token& token = tokens[v];
                    container.emplace_back(container.size(), token, index, left);
                    v = numeric_limits<size_t>::max();
                    left = LeftTokenOperator::NONE;
                }
                break;
            case '~':
                if(v == numeric_limits<size_t>::max()) {
                    // the reverse complement operator is specified before the token
                    left = LeftTokenOperator::REVERSE_COMPLEMENT;
                } else {
                    throw ConfigurationError(string("illegal right hand side operator in transform ") + c);
                }
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9': {
                if(v == numeric_limits<size_t>::max()) {
                    v = c - '0';
                } else {
                    v = v * 10 + (c - '0');
                }
                break;
            };
            default:
                throw ConfigurationError(string("illegal character in transform ") + c);
        }
        if(c == '\0') { break; }
        position++;
    }
};
void Environment::load_barcode_tolerance() {
    // multiplex barcode tolerance should have exactly the same number of elements
    // as the total number of multiplex barcode segments
    vector<uint8_t> provided(multiplex_barcode_tolerance);
    multiplex_barcode_tolerance.clear();
    for(size_t i = 0; i < total_multiplex_barcode_segments; i++) {
        if(i < provided.size()) {
            multiplex_barcode_tolerance.push_back(provided[i]);
        } else {
            multiplex_barcode_tolerance.push_back(numeric_limits<uint8_t>::max());
        }
    }
};
void Environment::load_thread_model() {
    if(transforms == numeric_limits<uint32_t>::max()) {
        transforms = threads;
    }
};
void Environment::load_undetermined() {
    /*  if there is only one channel and it has no barcode but its not marked undetermined
        we can logically mark it undetermined */
    if(channel_specifications.size() == 1) {
        ChannelSpecification* specification = channel_specifications[0];
        if(specification->multiplex_barcode.empty() && !specification->undetermined) {
            specification->undetermined = true;
        }
    }

    /*  if at most one channel is marked undetermined, set the undetermined pointer
        otherwise throw a validation exception */
    for(const auto specification : channel_specifications) {
        if(specification->undetermined) {
            if(undetermined == NULL) {
                undetermined = specification;
            } else {
                throw ConfigurationError("only one channel can accept undetermined reads");
            }
        }
    }
};
void Environment::load_multiplex_barcodes() {
    /*  compute a pairwise hamming distance metric for each barcode set
        and validate all barcodes obey the length constrains imposed by the transforms */
    for(size_t i = 0; i < total_multiplex_barcode_segments; i++) {
        Distance metric;
        size_t size = multiplex_barcode_length[i];
        for(const auto specification : channel_specifications) {
            if(!specification->multiplex_barcode.empty()) {
                if(specification->multiplex_barcode.size(i) == size) {
                    metric.add(specification->multiplex_barcode.iupac_ambiguity(i));
                } else {
                    string message("barcode ");
                    message.append(to_string(i));
                    message.append(" in channel ");
                    message.append(to_string(specification->index));
                    message.append(" is ");
                    message.append(to_string(specification->multiplex_barcode.size(i)));
                    message.append(" nucleotides long, expected ");
                    message.append(to_string(size));
                    throw ConfigurationError(message);
                }
            }
        }
        metric.load();
        multiplex_barcode_set_distance.emplace_back(metric);
    }

    /*  compute a pairwise hamming distance metric for the concatenated barcode */
    for(const auto specification : channel_specifications){
        if(!specification->multiplex_barcode.empty()) {
            multiplex_barcode_distance.add(specification->multiplex_barcode.iupac_ambiguity());
        }
    }
    multiplex_barcode_distance.load();

    /*  multiplex barcode tolerance can not be bigger than the number of correctable errors */
    for(size_t i = 0; i < total_multiplex_barcode_segments; i++) {
        multiplex_barcode_tolerance[i] = MIN(multiplex_barcode_tolerance[i], multiplex_barcode_set_distance[i].shannon_bound());
    }

    /*  construct a null barcode for the undetermined channel */
    if(undetermined != NULL && !multiplex_barcode_distance.empty()) {
        for(size_t i = 0; i < multiplex_barcode_set_distance.size(); i++) {
            Distance& metric = multiplex_barcode_set_distance[i];
            string pattern(metric.width(), '=');
            undetermined->multiplex_barcode.fill(i, pattern.c_str(), pattern.length());
        }
    }
};
void Environment::load_prior() {
    double total = 0;
    double count = 0;
    double specified = 0;
    for(const auto specification : channel_specifications) {
        if(!specification->undetermined) {
            count++;
            if(specification->concentration > 0) {
                total += specification->concentration;
                specified++;
            }
        }
    }
    if(specified > 0) {
        if(specified == count) {
            double factor = (1.0 - noise) / total;
            for(auto specification : channel_specifications) {
                specification->concentration *= factor;
            }
        } else {
            throw ConfigurationError("inconsistent channel concentration");
        }
    } else {
        /* if no concentrations were given assume a uniform distribution */
        double factor = (1.0 - noise) / count;
        for(auto specification : channel_specifications) {
            specification->concentration = factor;
        }
    }
    random_word_probability =  1.0 / double(pow(4, (multiplex_barcode_width())));;
    adjusted_noise_probability  = noise * random_word_probability;;
};
void Environment::load_channels() {
    HeadRGAtom default_rg;
    if(platform != Platform::UNKNOWN)   default_rg.set_platform(platform);
    if(!facility.empty())               kputsn(facility.c_str(), facility.size(), &default_rg.CN);
    if(!platform_model.empty())         kputsn(platform_model.c_str(), platform_model.size(), &default_rg.PM);
    if(!production_date.empty())        kputsn(production_date.c_str(), production_date.size(), &default_rg.DT);
    if(!insert_size.empty())            kputsn(insert_size.c_str(), insert_size.size(), &default_rg.PI);

    for(auto specifications : channel_specifications) {
        if(specifications->rg.ID.l > 0) {
            auto record = read_group_by_id.find(specifications->rg);
            if (record != read_group_by_id.end()) {
                specifications->rg.expand(*record->second);
            }
            specifications->rg.expand(default_rg);
            specifications->TC = total_output_segments;
            specifications->disable_quality_control = disable_quality_control;
            specifications->long_read = long_read;
            specifications->include_filtered = include_filtered;
            specifications->multiplex_barcode.set_threshold(masking_threshold);
            specifications->multiplex_barcode.set_tolerance(multiplex_barcode_tolerance);
        }
    }
};
void Environment::load_input_specification() {
    unordered_map< URL, size_t > referenced_input;

    for(const auto& url : input_urls) {
        referenced_input[url]++;
    }

    for(const auto& record : referenced_input) {
        const URL& url = record.first;
        FeedSpecification* feed = discover_feed(url, IoDirection::IN);
        feed->set_resolution(record.second);
        feed->set_capacity(buffer_capacity * feed->resolution);
    }
};
void Environment::load_output_specification() {
    unordered_map< URL, unordered_map< ChannelSpecification*, size_t > > referenced_output;

    for(const auto specifications : channel_specifications) {
        for(auto& url : specifications->output_urls) {
            referenced_output[url][specifications]++;
        }
    }

    for(const auto& url_record : referenced_output) {
        const URL& url = url_record.first;
        switch(url.type()) {
            case FormatType::SAM:
            case FormatType::BAM:
            case FormatType::CRAM:
            case FormatType::FASTQ: {
                size_t resolution = 0;
                for(const auto& channel_record : url_record.second) {
                    if(resolution == 0) {
                        resolution = channel_record.second;
                    } else if(resolution != channel_record.second) {
                        throw ConfigurationError("inconsistent resolution for " + string(url));
                    }
                }

                FeedSpecification* feed = discover_feed(url, IoDirection::OUT);
                feed->set_resolution(resolution);
                feed->set_capacity(buffer_capacity * feed->resolution);
                feed->register_pg(pg);

                for(const auto& channel_record : url_record.second) {
                    ChannelSpecification* channel = channel_record.first;
                    feed->register_rg(channel->rg);
                }
                break;
            };
            default:
                throw ConfigurationError("unknown format for " + string(url));
                break;
        }
    }
};
void Environment::validate_io() {
    /*  verify no URL is used for both input and output */
    for(const auto& output_record : output_feed_specification_by_url) {
        const URL& url = output_record.first;
        auto input_record = input_feed_specification_by_url.find(url);
        if (input_record != input_feed_specification_by_url.end()) {
            throw ConfigurationError("URL " + string(url) + " is used for both input and output");
        }
    }
};
void Environment::probe() {
    for(const auto& record : input_feed_specification_by_url) {
        FeedSpecification* feed = record.second;
        feed->probe();
    }
    for(const auto& record : output_feed_specification_by_url) {
        FeedSpecification* feed = record.second;
        feed->probe();
    }
};
void Environment::calibrate(const URL& url) {
    for(size_t i = 0; i < total_input_segments; i++) {
        input_urls.push_back(url);
        token_patterns.emplace_back(to_string(i) + "::");
        template_patterns.emplace_back(to_string(i));
    }
};
void Environment::validate() {
    /* validate phred range */
    if(input_phred_offset > MAX_PHRED_VALUE || input_phred_offset < MIN_PHRED_VALUE) {
        throw ConfigurationError("input phred offset out of range " + to_string(input_phred_offset));
    }
    if(output_phred_offset > MAX_PHRED_VALUE || output_phred_offset < MIN_PHRED_VALUE) {
        throw ConfigurationError("output phred offset out of range " + to_string(output_phred_offset));
    }

    switch (action) {
        case ProgramAction::DEMULTIPLEX: {

            /* leading_segment_index must reference an input segment */
            if (leading_segment_index >= total_input_segments) {
                throw ConfigurationError("invalid leading segment index " + to_string(leading_segment_index));
            }
            if (confidence < 0 || confidence > 1) {
                throw ConfigurationError("confidence value " + to_string(confidence) + " not between 0 and 1");
            }
            if (noise < 0 || noise > 1) {
                throw ConfigurationError("noise value " + to_string(noise) + " not between 0 and 1");
            }
            break;
        };
        case ProgramAction::QUALITY: {
            break;
        };
        default: break;
    }
};
FeedSpecification* Environment::discover_feed(const URL& url, const IoDirection& direction) {
    FeedSpecification* specification = NULL;
    switch(direction) {
        case IoDirection::IN: {
            auto record = input_feed_specification_by_url.find(url);
            if (record != input_feed_specification_by_url.end()) {
                specification = record->second;
            } else {
                specification = new FeedSpecification (
                    direction,
                    input_feed_specification_by_url.size(),
                    url,
                    platform,
                    input_phred_offset);
                input_feed_specification_by_url.emplace(make_pair(URL(url), specification));
            }
            break;
        };
        case IoDirection::OUT: {
            auto record = output_feed_specification_by_url.find(url);
            if (record != output_feed_specification_by_url.end()) {
                specification = record->second;
            } else {
                specification = new FeedSpecification (
                    direction,
                    output_feed_specification_by_url.size(),
                    url,
                    platform,
                    output_phred_offset);
                output_feed_specification_by_url.emplace(make_pair(URL(url), specification));
            }
            break;
        };
    }
    return specification;
};
