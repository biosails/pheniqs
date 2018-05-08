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

#ifndef PHENIQS_INTERFACE_H
#define PHENIQS_INTERFACE_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>

#include "error.h"
#include "url.h"
#include "json.h"

using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::list;
using std::size_t;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::unordered_map;
using std::invalid_argument;
using std::out_of_range;
using std::min;
using std::max;
using std::abs;
using std::istreambuf_iterator;

enum class ProgramAction : uint8_t {
    UNKNOWN,
    DEMULTIPLEX,
    QUALITY,
};
void to_string(const ProgramAction& value, string& result);
bool from_string(const char* value, ProgramAction& result);
bool from_string(const string& value, ProgramAction& result);
ostream& operator<<(ostream& o, const ProgramAction& value);
void encode_key_value(const string& key, const ProgramAction& value, Value& container, Document& document);

enum class ParameterType : uint8_t {
    BOOLEAN,
    INTEGER,
    DECIMAL,
    STRING,
    URL,
    DIRECTORY,
    UNKNOWN
};
void to_string(const ParameterType& value, string& result);
bool from_string(const char* value, ParameterType& result);
bool from_string(const string& value, ParameterType& result);
ostream& operator<<(ostream& o, const ParameterType& value);
void encode_key_value(const string& key, const ParameterType& value, Value& container, Document& document);

class Layout {
public:
    const size_t max_line_width;
    const int option_indent;
    int max_action_name;
    const int option_handle_spacing;
    const int option_choice_indent;
    Layout() :
        max_line_width(80),
        option_indent(2),
        max_action_name(0),
        option_handle_spacing(4),
        option_choice_indent(0) {
    };
    int complement_handle(const int& max_option_handle, const int& handle_length) const {
        return max_option_handle - handle_length + option_handle_spacing;
    };
    int complement_action(const int& action_name_length) const {
        return max_action_name - action_name_length + option_handle_spacing;
    };
    int indent_choice(const int& max_option_handle) const {
        return max_option_handle + option_handle_spacing + option_indent + option_choice_indent;
    };
};
class Prototype {
public:
    ParameterType type;
    string name;
    string help;
    string meta;
    uint64_t cardinality;
    bool plural;
    bool mandatory;
    bool positional;
    vector< string > handles;
    vector< string > choices;
    Prototype(const Value& node);
    Prototype();
    ostream& print_help(ostream& o, const int& max_option_handle, const Layout& layout) const;
    int handle_length() const;
    bool is_choice() const;
};
class Argument {
friend class CommandLine;
friend ostream& operator<<(ostream& o, const Argument& argument);

public:
    Argument(const Prototype* prototype);
    ~Argument();
    void set_value(const bool value);
    void set_value(const int64_t value);
    void set_value(const double value);
    void set_value(const char* value);
    void set_value(const URL value);
    bool* get_boolean() const;
    int64_t* get_integer() const;
    double* get_decimal() const;
    string* get_string() const;
    URL* get_url() const;
    URL* get_directory() const;
    list< URL >* get_directory_array() const;
    list< int64_t >* get_integer_array() const;
    list< double >* get_decimal_array() const;
    list< string >* get_string_array() const;
    list< URL >* get_url_array() const;
    uint64_t cardinality() const;
    bool satisfied() const;

private:
    bool assigned;
    const Prototype prototype;
    union {
        bool* boolean_value;
        int64_t* integer_value;
        double* decimal_value;
        string* string_value;
        URL* url_value;
        list< int64_t >* integer_array_value;
        list< double >* decimal_array_value;
        list< string >* string_array_value;
        list< URL >* url_array_value;
    };
};
class Action {
friend class CommandLine;

public:
    string name;
    string description;
    string epilog;
    Action(const Value& node, bool root=false);
    ~Action();
    inline int name_length() const {
        return static_cast< int >(name.size());
    };
    const Value* default_value();
    ostream& print_usage(ostream& o, const string& application_name, const Layout& layout) const;
    ostream& print_description_element(ostream& o, const Layout& layout) const;
    ostream& print_epilog_element(ostream& o, const Layout& layout) const;
    ostream& print_help(ostream& o, const string& application_name, const Layout& layout) const;
    bool is_root();

private:
    const bool root;
    int max_option_handle;
    const Value* default_value_node;
    vector< Prototype* > optional_order;
    vector< Prototype* > positional_order;
    unordered_map< string, Prototype* > option_handle_lookup;
    unordered_map< string, Prototype* > option_name_lookup;
};
class CommandLine {
public:
    const size_t argc;
    const char** argv;
    const string application_name;
    const string application_version;
    const string full_command;
    const URL working_directory;

    CommandLine(const int argc, const char** argv);
    ~CommandLine();
    const Document& instruction() const;
    const string& name() const;
    Argument* get(const string& name) const;
    ProgramAction get_selected_action() const;
    bool has_argument(const string& name);
    bool* get_boolean(const string& name) const;
    int64_t* get_integer(const string& name) const;
    double* get_decimal(const string& name) const;
    string* get_string(const string& name) const;
    URL* get_url(const string& name) const;
    URL* get_directory(const string& name) const;
    list< int64_t >* get_integer_plural(const string& name) const;
    list< double >* get_decimal_plural(const string& name) const;
    list< string >* get_string_plural(const string& name) const;
    list< URL >* get_url_plural(const string& name) const;
    list< URL >* get_directory_plural(const string& name) const;
    bool help_triggered() const;
    ostream& print_help(ostream& o) const;
    bool version_triggered() const;
    ostream& print_version(ostream& o) const;
    const Value* default_value() const;

private:
    Document _configuration;
    Document _instruction;
    Action* command;
    Action* selected;
    Layout layout;
    vector< Action* > action_order;
    unordered_map< string, Action* > action_name_lookup;
    vector< Argument* > argument_order;
    unordered_map< string, Argument* > argument_name_lookup;

    Argument* parse_argument(const Prototype* prototype, size_t& index);
    Argument* decode_optional(Action* action, size_t& index, const string& handle, bool& positional, bool composite=false);
    Argument* decode_positional(Action* action, size_t& index, const size_t& position);
    void decode();
    void validate();
    void load_action(size_t& position);
    Argument* get_argument(const Prototype* prototype);
    ostream& print_version_element(ostream& o, const Layout& layout) const;
    ostream& print_action_element(ostream& o, const Layout& layout) const;
    Document* load_default_instruction();
    Document* load_instruction_from_configuration_file();
    Document* load_instruction_from_command_line(const Value& base);
    void load_instruction();
};
#endif /* PHENIQS_INTERFACE_H */
