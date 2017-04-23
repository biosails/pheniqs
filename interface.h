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

#ifndef PHENIQS_INTERFACE_H
#define PHENIQS_INTERFACE_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <vector>
#include <unordered_map>

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/error/en.h>

#include "constant.h"
#include "error.h"

using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::unordered_map;
using std::invalid_argument;
using std::out_of_range;

using rapidjson::Document;
using rapidjson::Value;
using rapidjson::SizeType;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

class Dimension {
    public:
        const size_t max_line_width;
        const size_t option_indent;
        size_t max_action_name;
        const size_t option_handle_spacing;
        const size_t option_choice_indent;
        Dimension() :
            max_line_width(80),
            option_indent(2),
            max_action_name(0),
            option_handle_spacing(4),
            option_choice_indent(0) {
        };
        size_t complement_handle(const size_t& max_option_handle, const size_t& handle_length) const {
            return max_option_handle - handle_length + option_handle_spacing;
        };
        size_t complement_action(const size_t& action_name_length) const {
            return max_action_name - action_name_length + option_handle_spacing;
        };
        size_t indent_choice(const size_t& max_option_handle) const {
            return max_option_handle + option_handle_spacing + option_indent + option_choice_indent;
        };
};
class Prototype {
    public:
        ParameterType type;
        string name;
        string help;
        string meta;
        size_t cardinality;
        bool plural;
        bool mandatory;
        bool positional;
        vector<string> handles;
        vector<string> choices;
        Prototype() :
            cardinality(0),
            plural(false),
            mandatory(false),
            positional(false) {
        };
        ostream& print_help(ostream& o, const size_t& max_option_handle, const Dimension& dimension) const {
            o << setw(dimension.option_indent) << ' ';
            if(!positional) {
                for(size_t i = 0; i < handles.size(); i++) {
                    if(i > 0) {
                        o << ", ";
                    }
                    const string& handle = handles[i];
                    if(handle.length() == 1) {
                        o << '-' << handle;
                    } else {
                        o << "--" << handle;
                    }
                }
                if(meta.length() > 0) {
                    o << ' ';
                }
            }
            o << meta;
            o << setw(dimension.complement_handle(max_option_handle, handle_length())) << ' ';
            o << help << endl;
            return o;
        };
        size_t handle_length() const {
            size_t length = 0;
            if(!positional) {
                for(const auto& handle : handles) {
                    length += handle.length();
                    length += (handle.length() > 1 ? 4 : 3);
                }
                length--;
            }
            if(meta.length() > 0) {
                length++;
                length += meta.length();
            }
            return length;
        };
        bool is_choice() const {
            return choices.size() > 0;
        };
};
class Argument {
    friend class CommandLine;
    friend ostream& operator<<(ostream& o, const Argument& argument);

    public:
        Argument(const Prototype* prototype) :
            assigned(false),
            prototype(prototype) {
            if(prototype->plural) {
                switch(prototype->type) {
                    case ParameterType::boolean:
                        // does not make sense
                        break;
                    case ParameterType::integer:
                        integer_array_value = new vector<long>();
                        break;
                    case ParameterType::decimal:
                        decimal_array_value = new vector<double>();
                        break;
                    case ParameterType::string:
                        string_array_value = new vector<string>();
                        break;
                }
            } else {
                switch(prototype->type) {
                    case ParameterType::boolean:
                        boolean_value = new bool(false);
                        break;
                    case ParameterType::integer:
                        integer_value = new long(0);
                        break;
                    case ParameterType::decimal:
                        decimal_value = new double(0);
                        break;
                    case ParameterType::string:
                        string_value = new string();
                        break;
                }
            }
        };
        ~Argument() {
            if(prototype->plural) {
                switch(prototype->type) {
                    case ParameterType::boolean:
                        // does not make sense
                        break;
                    case ParameterType::integer:
                        delete integer_array_value;
                        break;
                    case ParameterType::decimal:
                        delete decimal_array_value;
                        break;
                    case ParameterType::string:
                        delete string_array_value;
                        break;
                }
            } else {
                switch(prototype->type) {
                    case ParameterType::boolean:
                        delete boolean_value;
                        break;
                    case ParameterType::integer:
                        delete integer_value;
                        break;
                    case ParameterType::decimal:
                        delete decimal_value;
                        break;
                    case ParameterType::string:
                        delete string_value;
                        break;
                }
            }
        };
        void set_value(const bool value) {
            assigned = true;
            *boolean_value = value;
        };
        void set_value(const long value) {
            assigned = true;
            if(!prototype->plural) {
                *integer_value = value;
            } else {
                integer_array_value->push_back(value);
            }
        };
        void set_value(const double value) {
            assigned = true;
            if(!prototype->plural) {
                *decimal_value = value;
            } else {
                decimal_array_value->push_back(value);
            }
        };
        void set_value(const char* value) {
            assigned = true;
            if(!prototype->plural) {
                string_value->assign(value);
            } else {
                string_array_value->emplace_back(value);
            }
        };
        bool* get_boolean() const {
            bool* value = NULL;
            if(prototype->type == ParameterType::boolean) {
                value = boolean_value;
            }
            return value;
        };
        long* get_integer() const {
            long* value = NULL;
            if(prototype->type == ParameterType::integer) {
                value = integer_value;
            }
            return value;
        };
        double* get_decimal() const {
            double* value = NULL;
            if(prototype->type == ParameterType::decimal) {
                value = decimal_value;
            }
            return value;
        };
        string* get_string() const {
            string* value = NULL;
            if(prototype->type == ParameterType::string) {
                value = string_value;
            }
            return value;
        };
        vector<long>* get_integer_array() const {
            vector<long>* value = NULL;
            if(prototype->type == ParameterType::integer && prototype->plural) {
                value = integer_array_value;
            }
            return value;
        };
        vector<double>* get_decimal_array() const {
            vector<double>* value = NULL;
            if(prototype->type == ParameterType::decimal && prototype->plural) {
                value = decimal_array_value;
            }
            return value;
        };
        vector<string>* get_string_array() const {
            vector<string>* value = NULL;
            if(prototype->type == ParameterType::string && prototype->plural) {
                value = string_array_value;
            }
            return value;
        };
        size_t cardinality() const {
            size_t value = 0;
            if(!prototype->plural) {
                value = assigned ? 1 : 0;
            } else {
                switch(prototype->type) {
                    case ParameterType::boolean:
                        // does not make sense
                        break;
                    case ParameterType::integer:
                        value = integer_array_value->size();
                        break;
                    case ParameterType::decimal:
                        value = decimal_array_value->size();
                        break;
                    case ParameterType::string:
                        value = string_array_value->size();
                        break;
                }
            }
            return value;
        };
        bool satisfied() const {
            bool value = true;
            if(prototype->plural) {
                if(prototype->cardinality == 0 || prototype->cardinality > cardinality()) {
                    value = false;
                }
            } else {
                if(!assigned) {
                    value = false;
                }
            }
            return value;
        };

    private:
        bool assigned;
        const Prototype* prototype;
        union {
            bool* boolean_value;
            long* integer_value;
            double* decimal_value;
            string* string_value;
            vector<long>* integer_array_value;
            vector<double>* decimal_array_value;
            vector<string>* string_array_value;
        };
};
class Action {
    friend class CommandLine;

    public:
        string name;
        string description;
        string epilog;
        Action(bool root=false) :
            root(root),
            max_option_handle(0) {
        };
        ~Action() {
            for(auto prototype : optional_order) {
                delete prototype;
            }
            for(auto prototype : positional_order) {
                delete prototype;
            }
        };
        bool is_root() {
            return root;
        };
        void append_option(Prototype* prototype) {
            if(!prototype->name.empty()) {
                if(option_name_lookup.count(prototype->name) == 0) {
                    option_name_lookup[prototype->name] = prototype;
                    for(const auto& handle : prototype->handles) {
                        if(option_handle_lookup.count(handle) == 0) {
                            option_handle_lookup[handle] = prototype;
                        } else {
                            throw ConfigurationError("redefining handle " + handle + " in option " + prototype->name);
                        }
                    }

                    // options without a handle are positional
                    if(prototype->handles.empty()) {
                        prototype->positional = true;
                    }

                    // default handle meta
                    if(!prototype->meta.length()) {
                        switch(prototype->type) {
                            case ParameterType::boolean:
                                prototype->meta.clear();
                                break;
                            case ParameterType::integer:
                                prototype->meta.assign("INT");
                                break;
                            case ParameterType::decimal:
                                prototype->meta.assign("FLOAT");
                                break;
                            case ParameterType::string:
                                prototype->meta.assign("STRING");
                                break;
                        }
                    }

                    if(prototype->positional) {
                        positional_order.push_back(prototype);
                    } else {
                        optional_order.push_back(prototype);
                    }
                    max_option_handle = MAX(max_option_handle, prototype->handle_length());
                } else {
                    throw ConfigurationError("redefining " + prototype->name);
                }
            } else {
                throw ConfigurationError("missing a name");
            }
        };
        ostream& print_usage(ostream& o, const string& application_name, const Dimension& dimension) const {
            string buffer;
            string block;
            buffer.append("Usage : ");
            buffer.append(application_name);
            if(!root) {
                buffer.append(" ");
                buffer.append(name);
            }
            size_t indent = buffer.length();
            for(const auto prototype : optional_order) {
                block.append(" ");
                if(!prototype->mandatory) {
                    block.append("[");
                }
                if(prototype->handles[0].length() == 1) {
                    block.append("-");
                } else {
                    block.append("--");
                }
                switch(prototype->type) {
                    case ParameterType::boolean:
                        block.append(prototype->handles[0]);
                        break;
                    case ParameterType::integer:
                        block.append(prototype->handles[0]);
                        block.append(" ");
                        block.append(prototype->meta);
                        break;
                    case ParameterType::decimal:
                        block.append(prototype->handles[0]);
                        block.append(" ");
                        block.append(prototype->meta);
                        break;
                    case ParameterType::string:
                        block.append(prototype->handles[0]);
                        block.append(" ");
                        if(prototype->is_choice()) {
                            for(size_t i = 0; i < prototype->choices.size(); i++) {
                                if(i > 0) {
                                    block.append("|");
                                }
                                block.append(prototype->choices[i]);
                            }
                        } else {
                            block.append(prototype->meta);
                        }
                        break;
                }
                if(!prototype->mandatory) {
                    block.append("]");
                }
                if(buffer.length() + block.length() > dimension.max_line_width) {
                    o << buffer << endl << setw(indent) << ' ';
                    buffer.clear();
                }
                buffer.append(block);
                block.clear();
            }
            for(const auto prototype : positional_order) {
                block.append(" ");
                if(!prototype->mandatory) {
                    block.append("[");
                }
                switch(prototype->type) {
                    case ParameterType::boolean:
                        break;
                    case ParameterType::integer:
                    case ParameterType::decimal:
                    case ParameterType::string: {
                        if(prototype->is_choice()) {
                            for(size_t i = 0; i < prototype->choices.size(); i++) {
                                if(i > 0) { block.append("|"); }
                                block.append(prototype->choices[i]);
                            }
                        } else {
                            block.append(prototype->meta);
                        }
                        break;
                    };
                }
                if(!prototype->mandatory) {
                    block.append("]");
                }
                if(prototype->plural) {
                    if(prototype->cardinality > 0) {
                        block.append("{");
                        block.append(to_string(prototype->cardinality));
                        block.append("}");
                    } else {
                        block.append("*");
                    }
                }
                if(buffer.length() + block.length() > dimension.max_line_width) {
                    o << buffer << endl << setw(indent) << ' ';
                    buffer.clear();
                }
                buffer.append(block);
                block.clear();
            }
            if(root) {
                buffer.append(" ACTION ...");
                if(buffer.length() + block.length() > dimension.max_line_width) {
                    o << buffer << endl << setw(indent) << ' ';
                    buffer.clear();
                }
                buffer.append(block);
                block.clear();
            }
            o << buffer;
            return o;
        };
        ostream& print_description_element(ostream& o, const Dimension& dimension) const {
            if(description.length() > 0) {
                o << description << endl << endl;
            }
            return o;
        };
        ostream& print_epilog_element(ostream& o, const Dimension& dimension) const {
            if(epilog.length() > 0) {
                o << endl << epilog << endl;
            }
            return o;
        };
        ostream& print_help(ostream& o, const string& application_name, const Dimension& dimension) const {
            print_usage(o, application_name, dimension);
            if(!positional_order.empty()) {
                o << endl << endl << "Positional:" << endl;
                for(const auto option : positional_order) {
                    option->print_help(o, max_option_handle, dimension);
                }
            }
            if(!optional_order.empty()) {
                o << endl << endl << "Optional:" << endl;
                for(const auto option : optional_order) {
                    option->print_help(o, max_option_handle, dimension);
                }
            }
            return o;
        };

    private:
        const bool root;
        size_t max_option_handle;
        vector< Prototype* > optional_order;
        vector< Prototype* > positional_order;
        unordered_map< string, Prototype* > option_handle_lookup;
        unordered_map< string, Prototype* > option_name_lookup;
};
class CommandLine {
    public:
        CommandLine(const char* content, const string& version) :
            argc(0),
            argv(NULL),
            version(version),
            command(NULL),
            selected(NULL) {
            load_configuration(content);
        };
        ~CommandLine() {
            for(auto action : action_order) {
                delete action;
            }
            for(auto argument : argument_order) {
                delete argument;
            }
        };
        void load(int argc, char** argv) {
            this->argc = size_t(argc);
            this->argv = argv;
            decode();
            validate();
        };
        string& name() {
            if(!command->name.empty()) {
                return command->name;
            } else {
                return application_name;
            }
        };
        string& get_full_command() {
            return full_command;
        };
        Argument* get(const string& name) {
            auto record = argument_name_lookup.find(name);
            if (record != argument_name_lookup.end()) {
                return record->second;
            } else {
                return NULL;
            }
        };
        string& get_selected_action() {
            return selected->name;
        };
        bool has_argument(const string& name) {
            return argument_name_lookup.count(name) > 0;
        };
        bool* get_boolean(const string& name) {
            bool* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_boolean();
            }
            return value;
        };
        long* get_integer(const string& name) {
            long* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_integer();
            }
            return value;
        };
        double* get_decimal(const string& name) {
            double* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_decimal();
            }
            return value;
        };
        string* get_string(const string& name) {
            string* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_string();
            }
            return value;
        };
        vector<long>* get_integer_plural(const string& name) {
            vector<long>* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_integer_array();
            }
            return value;
        };
        vector<double>* get_decimal_plural(const string& name) {
            vector<double>* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_decimal_array();
            }
            return value;
        };
        vector<string>* get_string_plural(const string& name) {
            vector<string>* value = NULL;
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_string_array();
            }
            return value;
        };
        bool help_triggered() {
            if(argc < 2) {
                return true;
            } else {
                bool* value = get_boolean("help");
                if(value != NULL && value) {
                    return true;
                } else {
                    return false;
                }
            }
        };
        ostream& print_help(ostream& o) {
            o << endl;
            print_version_element(o, dimension);
            if(!selected->root) {
                command->print_description_element(o, dimension);
            }
            selected->print_description_element(o, dimension);
            selected->print_help(o, name(), dimension);
            if(selected->root) {
                print_action_element(o, dimension);
            }
            selected->print_epilog_element(o, dimension);
            o << endl;
            return o;
        };
        bool version_triggered() {
            bool* value = get_boolean("version");
            if(value != NULL && value) {
                return true;
            } else {
                return false;
            }
        };
        ostream& print_version(ostream& o) {
            print_version_element(o, dimension);
            return o;
        };
    private:
        size_t argc;
        char** argv;
        string full_command;
        string version;
        Action* command;
        Action* selected;
        string application_name;
        Dimension dimension;
        vector< Action* > action_order;
        unordered_map< string, Action* > action_name_lookup;
        vector< Argument* > argument_order;
        unordered_map< string, Argument* > argument_name_lookup;
        void load_configuration(const char* content) {
            Document document;
            if (!document.Parse(content).HasParseError()) {
                if (document.IsObject()) {
                    command = load_action_node(document, true);
                    Value::ConstMemberIterator element = document.FindMember("action");
                    if (element != document.MemberEnd()) {
                        if(element->value.IsArray()) {
                            for (SizeType i = 0; i < element->value.Size(); i++) {
                                try {
                                    add_action(load_action_node(element->value[i], false));
                                } catch(ConfigurationError& e) {
                                    e.message += " in action " + to_string(i);
                                    throw e;
                                }
                            }
                        } else { throw ConfigurationError("incorrect action syntax"); }
                    }
                } else { throw ConfigurationError("incorrect action syntax"); }
            } else {
                throw ConfigurationError(string(GetParseError_En(document.GetParseError())) + " at position " + to_string(document.GetErrorOffset()));
            }
        };
        Action* load_action_node(const Value& node, bool root) {
            Action* action = NULL;
            if (node.IsObject()) {
                action = new Action(root);
                Value::ConstMemberIterator element;
                element = node.FindMember("name");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        action->name.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect syntax"); }
                }
                element = node.FindMember("description");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        action->description.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect syntax"); }
                }
                element = node.FindMember("epilog");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        action->epilog.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect syntax"); }
                }
                element = node.FindMember("option");
                if (element != node.MemberEnd()) {
                    if(element->value.IsArray()) {
                        for (SizeType i = 0; i < element->value.Size(); i++) {
                            action->append_option(load_prototype_node(element->value[i]));
                        }
                    } else { throw ConfigurationError("incorrect syntax"); }
                }
            } else {
                throw ConfigurationError("incorrect action syntax");
            }
            return action;
        };
        Prototype* load_prototype_node(const Value& node) {
            Prototype* prototype = NULL;
            if (node.IsObject()) {
                prototype = new Prototype();
                Value::ConstMemberIterator element;
                element = node.FindMember("name");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        prototype->name.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("help");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        prototype->help.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("meta");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        prototype->meta.assign(element->value.GetString(), element->value.GetStringLength());
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("cardinality");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        if(element->value.GetStringLength() == 1 && *(element->value.GetString()) == '*') {
                            prototype->plural = true;
                        } else { throw ConfigurationError("incorrect prototype syntax"); }
                    } else if(element->value.IsUint()) {
                        prototype->plural = true;
                        prototype->cardinality = element->value.GetUint64();
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("mandatory");
                if (element != node.MemberEnd()) {
                    if(element->value.IsBool()) {
                        prototype->mandatory = element->value.GetBool();
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("type");
                if (element != node.MemberEnd()) {
                    if(element->value.IsString()) {
                        const Value::Ch* value = element->value.GetString();
                        if(!strcmp(value, "boolean")) {
                            prototype->type = ParameterType::boolean;
                        } else if (!strcmp(value, "integer")) {
                            prototype->type = ParameterType::integer;
                        } else if (!strcmp(value, "decimal")) {
                            prototype->type = ParameterType::decimal;
                        } else if (!strcmp(value, "string")) {
                            prototype->type = ParameterType::string;
                        } else {
                            throw ConfigurationError("unknown argument type " + string(value));
                        }
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("handle");
                if (element != node.MemberEnd()) {
                    if(element->value.IsArray()) {
                        for (SizeType i = 0; i < element->value.Size(); i++) {
                            if(element->value[i].IsString()) {
                                SizeType length = element->value[i].GetStringLength();
                                const Value::Ch* value = element->value[i].GetString();
                                if(length == 2 && value[0] == '-' && value[1] != '-') {
                                    // short
                                    prototype->handles.emplace_back(value, 1, length - 1);
                                } else if(length > 3 && value[0] == '-' && value[1] == '-' && value[2] != '-') {
                                    // long
                                    prototype->handles.emplace_back(value, 2, length - 2);
                                } else { throw ConfigurationError("incorrect option handle syntax " + string(value, length)); }
                            } else { throw ConfigurationError("incorrect prototype syntax"); }
                        }
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
                element = node.FindMember("choice");
                if (element != node.MemberEnd()) {
                    if(element->value.IsArray()) {
                        for (SizeType i = 0; i < element->value.Size(); i++) {
                            if(element->value[i].IsString()) {
                                prototype->choices.emplace_back(element->value[i].GetString(), element->value[i].GetStringLength());
                            } else { throw ConfigurationError("incorrect prototype syntax"); }
                        }
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
            } else { throw ConfigurationError("incorrect prototype syntax"); }
            return prototype;
        };
        inline void add_action(Action* action) {
            if(!action->name.empty()) {
                if(!action_name_lookup.count(action->name)) {
                    action_order.push_back(action);
                    action_name_lookup[action->name] = action;
                    dimension.max_action_name = MAX(dimension.max_action_name, action->name.length());
                } else { throw ConfigurationError("redefining " + action->name); }
            } else { throw ConfigurationError("action missing a name"); }
        };
        inline Argument* parse_argument(const Prototype* prototype, size_t& index) {
            Argument* argument = NULL;
            try {
                switch(prototype->type) {
                    case ParameterType::boolean: {
                        argument = get_argument(prototype);
                        argument->set_value(true);
                        break;
                    };
                    case ParameterType::integer: {
                        string raw(argv[index]);
                        argument = get_argument(prototype);
                        argument->set_value(stol(raw));
                        break;
                    };
                    case ParameterType::decimal: {
                        string raw(argv[index]);
                        argument = get_argument(prototype);
                        argument->set_value(stod(raw));
                        break;
                    };
                    case ParameterType::string: {
                        argument = get_argument(prototype);
                        argument->set_value(argv[index]);
                        break;
                    };
                }
            } catch(invalid_argument& e) {
                throw CommandLineError("invalid value for argument " + prototype->name);
            } catch(out_of_range& e) {
                throw CommandLineError("out of range value for argument " + prototype->name);
            }
            return argument;
        };
        inline Argument* decode_optional(Action* action, size_t& index, const string& handle, bool& positional, bool composite=false) {
            Argument* argument = NULL;
            auto record = action->option_handle_lookup.find(handle);
            if (record != action->option_handle_lookup.end()) {
                Prototype* prototype = record->second;
                // if the option is not a boolean flag
                if(prototype->type != ParameterType::boolean) {
                    if(composite) {
                        // if the handle was grouped and not the last in the group it must be a boolean option
                        throw CommandLineError("argument " + prototype->name + " requires a value");
                    } else  if(index + 1 < argc) {
                        // if there is no next token or the next token is a noop than the option is missing a value
                        index++;
                        if(!strcmp(argv[index], "--")) {
                            positional = true;
                            index++;
                            throw CommandLineError("argument " + prototype->name + " requires a value");
                        }
                    } else { throw CommandLineError("argument " + prototype->name + " requires a value"); }
                }
                argument = parse_argument(prototype, index);
            }
            return argument;
        };
        inline Argument* decode_positional(Action* action, size_t& index, const size_t& position) {
            Argument* argument = NULL;
            if(position < action->positional_order.size()) {
                Prototype* prototype = action->positional_order[position];
                argument = parse_argument(prototype, index);
            } else {
                throw CommandLineError("too many positional arguments");
            }
            return argument;
        };
        inline void decode() {
            application_name.assign(argv[0]);
            load_full_command();
            size_t index = 1;
            size_t position = 0;
            bool positional = false;
            load_action(index);
            while(index < argc) {
                if(!strcmp(argv[index], "--")) {
                    // break signal encountered
                    if(positional) {
                        // if parsing positional arguments skip to the next positional
                        position++;
                    } else {
                        // if parsing optional switch to decoding positional
                        positional = true;
                    }
                } else {
                    if(!positional) {
                        char* key = argv[index];
                        size_t length = strlen(key);
                        if (key[0] == '-' && length > 1) {
                            key++;
                            if(*key == '-') {
                                if(length > 3) {
                                    // long
                                    key++;
                                    string handle(key);
                                    if(!decode_optional(selected, index, handle, positional)) {
                                        throw CommandLineError("unknown argument --" + handle);
                                    }
                                } else {
                                    throw CommandLineError("unknown argument " + string(argv[index]));
                                }
                            } else {
                                // short, potentially clustered
                                while(*key != '\0') {
                                    string handle(1, *key);
                                    if(!decode_optional(selected, index, handle, positional, *(key + 1) != '\0')) {
                                        throw CommandLineError("unknown argument -" + handle);
                                    }
                                    key++;
                                }
                            }
                        } else {
                            // first positional
                            positional = true;
                        }
                    }
                    if(positional) {
                        if(position < selected->positional_order.size()) {
                            Argument* argument = decode_positional(selected, index, position);
                            if(argument->satisfied()) {
                                position++;
                            }
                        }
                    }
                }
                index++;
            }
        };
        inline void validate() {
            if(!help_triggered()) {
                for(const auto& record : selected->option_name_lookup) {
                    Argument* argument = NULL;
                    const string& name = record.first;
                    const Prototype* prototype = record.second;

                    auto argument_record = argument_name_lookup.find(name);
                    if (argument_record != argument_name_lookup.end()) {
                        argument = argument_record->second;
                    }
                    if(prototype->mandatory && (argument == NULL || !argument->assigned)) {
                        throw CommandLineError("missing mandatory argument " + prototype->name);
                    } else if(prototype->plural && prototype->cardinality > 0 && (argument == NULL || prototype->cardinality != argument->cardinality())) {
                        throw CommandLineError("argument " + prototype->name + " requires exactly " + to_string(prototype->cardinality) + " values");
                    } else if(argument != NULL && prototype->is_choice()) {
                        switch(prototype->type) {
                            case ParameterType::boolean: {
                                break;
                            };
                            case ParameterType::integer: {
                                break;
                            };
                            case ParameterType::decimal: {
                                break;
                            };
                            case ParameterType::string: {
                                bool found = false;
                                string* value = argument->get_string();
                                if(value != NULL) {
                                    for(const auto& choice : prototype->choices) {
                                        if(!choice.compare(*value)) {
                                            found = true;
                                            break;
                                        }
                                    }
                                }
                                if(!found){
                                    throw CommandLineError("invalid value " + *value + " for option " + prototype->name);
                                }
                                break;
                            };
                        };
                    }

                }
            }
        };
        inline void load_full_command() {
            for(size_t i = 0; i < argc; i++) {
                full_command.append(argv[i]);
                full_command.append(" ");
            }
        };
        inline void load_action(size_t& position) {
            if (argc > 1) {
                for (size_t i = 1; i < argc; i++) {
                    if (argv[i][0] != '-') {
                        string name(argv[i]);
                        position = i + 1;
                        auto record = action_name_lookup.find(name);
                        if (record != action_name_lookup.end()) {
                            selected = record->second;
                        } else {
                            throw CommandLineError("unknown action " +  name);
                        }
                        break;
                    }
                }
            }
            if(selected == NULL) {
                selected = command;
            }
        };
        Argument* get_argument(const Prototype* prototype) {
            Argument* argument = NULL;
            auto record = argument_name_lookup.find(prototype->name);
            if (record != argument_name_lookup.end()) {
                argument = record->second;
            } else {
                argument = new Argument(prototype);
                argument_name_lookup[prototype->name] = argument;
                argument_order.push_back(argument);
            }
            return argument;
        };
        ostream& print_version_element(ostream& o, const Dimension& dimension) {
            o << name() << " version " << version << endl;
            return o;
        };
        ostream& print_action_element(ostream& o, const Dimension& dimension) {
            o << endl << "available action" << endl;
            for(const auto action : action_order) {
                o << setw(dimension.option_indent) << ' ';
                o << action->name;
                o << setw(dimension.complement_action(action->name.length())) << ' ';
                o << action->description;
                o << endl;
            }
            return o;
        };
};
#endif /* PHENIQS_INTERFACE_H */