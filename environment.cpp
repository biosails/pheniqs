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

#include "environment.h"

/*  Environment */
Environment::Environment(const int argc, const char** argv) :
    instruction(kObjectType),
    interface(argc, argv),
    _program_action(interface.get_selected_action()),
    _help_only(interface.help_triggered()),
    _version_only(interface.version_triggered()),
    _validate_only(false),
    _lint_only(false),
    _display_distance(false) {

    instruction.CopyFrom(interface.instruction(), instruction.GetAllocator());
    if(!_help_only && !_version_only) {
        switch (_program_action) {
            case ProgramAction::DEMULTIPLEX: {
                decode_value_by_key< bool >("validate only", _validate_only, instruction);
                decode_value_by_key< bool >("lint only", _lint_only, instruction);
                decode_value_by_key< bool >("display distance", _display_distance, instruction);
                break;
            };

            case ProgramAction::QUALITY: {
                break;
            };

            default: break;
        }
    }
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
    instruction.Accept(writer);
    o << buffer.GetString() << endl;
};
void Environment::print_instruction_validation(ostream& o) const {
    print_global_instruction(o);
    print_input_instruction(o);
    print_template_instruction(o);
    print_multiplex_instruction(o);
    print_molecular_instruction(o);
    print_splitseq_instruction(o);
};
void Environment::print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const {
    Value::ConstMemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
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
void Environment::print_codec_instruction(const Value& value, const bool& plural, ostream& o) const {
    if(!value.IsNull()) {
        if(plural) {
            int32_t index;
            if(decode_value_by_key< int32_t >("index", index, value)) {
                o << "  Decoder No." << to_string(index) << endl << endl;
            }
        }

        Algorithm algorithm(decode_value_by_key< Algorithm >("algorithm", value));
        o << "    Decoding algorithm                   " << algorithm << endl;

        print_codec_template(value, o);
        if(_display_distance) {
            CodecDistanceMetric metric(value);
            metric.describe(cout);
        }

        Value::ConstMemberIterator reference = value.FindMember("codec");
        if(reference != value.MemberEnd()) {
            if(reference->value.IsObject()) {
                for(const auto& record : reference->value.GetObject()) {
                    const Value& element(record.value);
                    int32_t index(decode_value_by_key< int32_t >("index", element));
                    o << "    Channel No." << index << endl;

                    string buffer;
                    if(decode_value_by_key< string >("ID", buffer, element)) { o << "        ID : " << buffer << endl; }
                    if(decode_value_by_key< string >("PU", buffer, element)) { o << "        PU : " << buffer << endl; }
                    if(decode_value_by_key< string >("LB", buffer, element)) { o << "        LB : " << buffer << endl; }
                    if(decode_value_by_key< string >("SM", buffer, element)) { o << "        SM : " << buffer << endl; }
                    if(decode_value_by_key< string >("DS", buffer, element)) { o << "        DS : " << buffer << endl; }
                    if(decode_value_by_key< string >("DT", buffer, element)) { o << "        DT : " << buffer << endl; }
                    if(decode_value_by_key< string >("PL", buffer, element)) { o << "        PL : " << buffer << endl; }
                    if(decode_value_by_key< string >("PM", buffer, element)) { o << "        PM : " << buffer << endl; }
                    if(decode_value_by_key< string >("CN", buffer, element)) { o << "        CN : " << buffer << endl; }
                    if(decode_value_by_key< string >("FO", buffer, element)) { o << "        FO : " << buffer << endl; }
                    if(decode_value_by_key< string >("KS", buffer, element)) { o << "        KS : " << buffer << endl; }
                    if(decode_value_by_key< string >("PI", buffer, element)) { o << "        PI : " << buffer << endl; }
                    if(decode_value_by_key< string >("FS", buffer, element)) { o << "        FS : " << buffer << endl; }
                    if(decode_value_by_key< string >("CO", buffer, element)) { o << "        CO : " << buffer << endl; }
                    double concentration;
                    if(decode_value_by_key< double >("concentration", concentration, element)) {
                        o << "        Concentration : " << concentration << endl;
                    }

                    int32_t segment_index(0);
                    list< URL > output;
                    if(decode_value_by_key< list< URL > >("output", output, element)) {
                        for(auto& url : output) {
                            o << "        Segment No." + to_string(segment_index) + "  : " << url << endl;
                            ++segment_index;
                        }
                    }
                    o << endl;
                }
            }
        }
    }
};
void Environment::print_feed_instruction(const Value::Ch* key, ostream& o) const {
    Value::ConstMemberIterator reference = instruction.FindMember(key);
    if(reference != instruction.MemberEnd()) {
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
void Environment::print_global_instruction(ostream& o) const {
    o << fixed << setprecision(19);
    o << "Environment " << endl << endl;
    o << "    Version                                     " << interface.application_version << endl;

    URL base_input_url;
    decode_value_by_key< URL >("base input url", base_input_url, instruction);
    o << "    Base input URL                              " << base_input_url << endl;

    URL base_output_url;
    decode_value_by_key< URL >("base input url", base_output_url, instruction);
    o << "    Base output URL                             " << base_output_url << endl;

    Platform platform;
    decode_value_by_key< Platform >("platform", platform, instruction);
    o << "    Platform                                    " << platform << endl;

    bool disable_quality_control;
    decode_value_by_key< bool >("disable quality control", disable_quality_control, instruction);
    o << "    Quality tracking                            " << (disable_quality_control ? "disabled" : "enabled") << endl;

    bool include_filtered;
    decode_value_by_key< bool >("include filtered", include_filtered, instruction);
    o << "    Include non PF reads                        " << (include_filtered ? "enabled" : "disabled") << endl;

    uint8_t input_phred_offset;
    decode_value_by_key< uint8_t >("input phred offset", input_phred_offset, instruction);
    o << "    Input Phred offset                          " << to_string(input_phred_offset) << endl;

    uint8_t output_phred_offset;
    decode_value_by_key< uint8_t >("output phred offset", output_phred_offset, instruction);
    o << "    Output Phred offset                         " << to_string(output_phred_offset) << endl;

    int32_t leading_segment_index;
    decode_value_by_key< int32_t >("leading segment index", leading_segment_index, instruction);
    o << "    Leading template segment                    " << to_string(leading_segment_index) << endl;

    int32_t buffer_capacity;
    decode_value_by_key< int32_t >("buffer capacity", buffer_capacity, instruction);
    o << "    Feed buffer capacity                        " << to_string(buffer_capacity) << endl;

    int32_t threads;
    decode_value_by_key< int32_t >("threads", threads, instruction);
    o << "    Threads                                     " << to_string(threads) << endl;
    o << endl;
};
void Environment::print_input_instruction(ostream& o) const {
    o << "Input " << endl << endl;

    int32_t input_segment_cardinality;
    if(decode_value_by_key< int32_t >("input segment cardinality", input_segment_cardinality, instruction)) {
        o << "    Input segment cardinality                   " << to_string(input_segment_cardinality) << endl;
    }

    list< URL > input_url_array;
    if(decode_value_by_key< list< URL > >("input", input_url_array, instruction)) {
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
void Environment::print_template_instruction(ostream& o) const {
    o << "Template" << endl << endl;

    int32_t output_segment_cardinality;
    if(decode_value_by_key< int32_t >("output segment cardinality", output_segment_cardinality, instruction)) {
        o << "    Output segment cardinality                  " << to_string(output_segment_cardinality) << endl;
    }

    Rule template_rule(decode_value_by_key< Rule >("template", instruction));
    o << endl;
    for(auto& token : template_rule.token_array) {
        o << "    Token No." << token.index << endl;
        o << "        Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
        o << "        Pattern       " << string(token) << endl;
        o << "        Description   ";
        o << token.description() << endl;
        o << endl;
    }
    o << "    Transformation" << endl;
    for(const auto& transform : template_rule.transform_array) {
        o << "        " << transform.description() << endl;
    }
    o << endl;
};
void Environment::print_codec_template(const Value& value, ostream& o) const {
    int32_t segment_cardinality;
    if(decode_value_by_key< int32_t >("segment cardinality", segment_cardinality, value)) {
        o << "    Segment cardinality                  " << to_string(segment_cardinality) << endl;
    }

    int32_t nucleotide_cardinality;
    if(decode_value_by_key< int32_t >("nucleotide cardinality", nucleotide_cardinality, value)) {
        o << "    Nucleotide cardinality               " << to_string(nucleotide_cardinality) << endl;
    }

    vector< int32_t > barcode_length;
    if(decode_value_by_key< vector< int32_t > >("barcode length", barcode_length, value)) {
        o << "    Barcode segment length               ";
        for(const auto& v : barcode_length) {
            o << to_string(v) << " ";
        }
        o << endl;
    }

    Rule rule(decode_value_by_key< Rule >("template", value));
    o << endl;
    for(auto& token : rule.token_array) {
        o << "    Token No." << token.index << endl;
        o << "        Length        " << (token.constant() ? to_string(token.length()) : "variable") << endl;
        o << "        Pattern       " << string(token) << endl;
        o << "        Description   ";
        o << token.description() << endl;
        o << endl;
    }
    o << "    Transform" << endl;
    for(const auto& transform : rule.transform_array) {
        o << "        " << transform.description() << endl;
    }
    o << endl;
};
void Environment::print_multiplex_instruction(ostream& o) const {
    print_codec_group_instruction("multiplex", "Mutliplexing", o);
    print_feed_instruction("output feed", o);
};
void Environment::print_molecular_instruction(ostream& o) const {
    print_codec_group_instruction("molecular", "Unique Molecular Identifier", o);
};
void Environment::print_splitseq_instruction(ostream& o) const {
    print_codec_group_instruction("splitseq", "SplitSEQ", o);
};
void Environment::load_pheniqs_pg() {
    // ks_put_string(interface.name().c_str(), interface.name().size(), &_pheniqs_pg.ID);
    // ks_put_string(interface.name().c_str(), interface.name().size(), &_pheniqs_pg.PN);
    // ks_put_string(interface.application_version.c_str(), interface.application_version.size(), &_pheniqs_pg.VN);
    // ks_put_string(interface.full_command.c_str(), interface.full_command.size(), &_pheniqs_pg.CL);
};
