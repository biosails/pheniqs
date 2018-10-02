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

#include "include.h"
#include "version.h"
#include "url.h"

class Layout;
class Prototype;
class Argument;
class Action;
class Interface;

enum class ParameterType : uint8_t {
    BOOLEAN,
    INTEGER,
    DECIMAL,
    STRING,
    URL,
    UNKNOWN
};
string to_string(const ParameterType& value);
bool from_string(const char* value, ParameterType& result);
bool from_string(const string& value, ParameterType& result);
ostream& operator<<(ostream& o, const ParameterType& value);
void encode_key_value(const string& key, const ParameterType& value, Value& container, Document& document);
template<> bool decode_value_by_key< ParameterType >(const Value::Ch* key, ParameterType& value, const Value& container);
template<> ParameterType decode_value_by_key< ParameterType >(const Value::Ch* key, const Value& container);

class Layout {
    public:
        const size_t max_line_width;
        const int option_indent;
        int max_action_name;
        int max_option_handle;
        const int option_handle_spacing;
        const int option_choice_indent;
        Layout() :
            max_line_width(80),
            option_indent(2),
            max_action_name(0),
            option_handle_spacing(4),
            option_choice_indent(0) {
        };
        int pad_option_handle(const int& handle_length) const {
            return max_option_handle - handle_length + option_handle_spacing;
        };
        int pad_action_name(const int& action_name_length) const {
            return max_action_name - action_name_length + option_handle_spacing;
        };
        void load_action(const Action& action);
};

class Prototype {
    public:
        const ParameterType type;
        const string name;
        string help;
        string meta;
        uint32_t cardinality;
        bool plural;
        bool mandatory;
        bool positional;
        vector< string > handles;
        vector< string > choices;
        Prototype(const Value& ontology);
        inline bool is_choice() const {
            return choices.size() > 0;
        };
        int handle_length() const;
        ostream& print_help(ostream& o, const Layout& layout) const;
};

class Argument {
    friend ostream& operator<<(ostream& o, const Argument& argument);

    public:
        Argument(const Prototype* prototype);
        ~Argument();
        inline const bool empty() const {
            return !assigned;
        };
        inline const bool plural() const {
            return prototype.plural;
        };
        inline const ParameterType type() const {
            return prototype.type;
        };
        inline void set_value(const bool value) {
            assigned = true;
            *boolean_value = value;
        };
        inline void set_value(const int64_t value) {
            assigned = true;
            if(!prototype.plural) {
                *integer_value = value;
            } else {
                integer_array_value->push_back(value);
            }
        };
        inline void set_value(const double value) {
            assigned = true;
            if(!prototype.plural) {
                *decimal_value = value;
            } else {
                decimal_array_value->push_back(value);
            }
        };
        inline void set_value(const char* value) {
            assigned = true;
            if(!prototype.plural) {
                string_value->assign(value);
            } else {
                string_array_value->emplace_back(value);
            }
        };
        inline void set_value(const URL value) {
            assigned = true;
            if(!prototype.plural) {
                *url_value = value;
            } else {
                url_array_value->emplace_back(value);
            }
        };
        bool* get_boolean() const;
        int64_t* get_integer() const;
        double* get_decimal() const;
        string* get_string() const;
        URL* get_url() const;
        list< int64_t >* get_integer_array() const;
        list< double >* get_decimal_array() const;
        list< string >* get_string_array() const;
        list< URL >* get_url_array() const;
        uint32_t cardinality() const;
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
    friend class Layout;
    friend bool encode_key_value(const string& key, const Action& value, Value& container, Document& document);
    public:
        const Value& ontology;
        const string name;
        const string description;
        Action(const Value& ontology, bool root=false);
        virtual ~Action();
        Document operation();
        virtual void load(const size_t argc, const char** argv);
        inline const bool is_root() const {
            return root;
        };
        inline int name_length() const {
            return static_cast< int >(name.size());
        };
        bool help_triggered() const {
            bool* value(get_boolean("help only"));
            return (value == NULL) ? false : *value;
        };
        bool version_triggered() const {
            bool* value(get_boolean("version only"));
            return (value == NULL) ? false : *value;
        };
        ostream& print_usage(ostream& o, const string& application_name, const Layout& layout) const;
        ostream& print_description(ostream& o, const Layout& layout) const;
        ostream& print_epilog(ostream& o, const Layout& layout) const;
        ostream& print_prolog(ostream& o, const Layout& layout) const;
        ostream& print_license(ostream& o, const Layout& layout) const;
        ostream& print_positional(ostream& o, const Layout& layout) const;
        ostream& print_optional(ostream& o, const Layout& layout) const;
        ostream& print_help(ostream& o, const string& application_name, const Layout& layout) const;
        ostream& print_action_dictionary_item(ostream& o, const Layout& layout) const;

    private:
        const bool root;
        list< Prototype* > optional_by_index;
        vector< Prototype* > positional_by_index;
        unordered_map< string, Prototype* > option_by_handle;
        unordered_map< string, Prototype* > option_by_name;
        list< Argument* > argument_by_index;
        unordered_map< string, Argument* > argument_by_name;
        Argument* get(const string& name) const {
            auto record = argument_by_name.find(name);
            if(record != argument_by_name.end()) {
                return record->second;
            } else {
                return NULL;
            }
        };
        bool* get_boolean(const string& name) const {
            bool* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_boolean();
            }
            return value;
        };
        int64_t* get_integer(const string& name) const {
            int64_t* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_integer();
            }
            return value;
        };
        double* get_decimal(const string& name) const {
            double* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_decimal();
            }
            return value;
        };
        string* get_string(const string& name) const {
            string* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_string();
            }
            return value;
        };
        URL* get_url(const string& name) const {
            URL* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_url();
            }
            return value;
        };
        list< int64_t >* get_integer_plural(const string& name) const {
            list< int64_t >* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_integer_array();
            }
            return value;
        };
        list< double >* get_decimal_plural(const string& name) const {
            list< double >* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_decimal_array();
            }
            return value;
        };
        list< string >* get_string_plural(const string& name) const {
            list< string >* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_string_array();
            }
            return value;
        };
        list< URL >* get_url_plural(const string& name) const {
            list< URL >* value(NULL);
            Argument* argument = get(name);
            if(argument) {
                value = argument->get_url_array();
            }
            return value;
        };
        Argument* get_argument(const Prototype* prototype);
        Argument* parse_argument(const size_t argc, const char** argv, const Prototype* prototype, size_t& index);
        Argument* load_optional_argument(const size_t argc, const char** argv, size_t& index, const string& handle, bool& positional, bool composite=false);
        Argument* load_positional_argument(const size_t argc, const char** argv, size_t& index, const size_t& position);
        Document load_document(const URL& url);
        void set_help_triggered();

    protected:
        virtual void parse_command_line(const size_t argc, const char** argv);
        virtual void validate_command_line();
};

class Interface {
    public:
        const size_t argc;
        const char** argv;
        const string application_name;
        const string application_version;
        const string full_command;
        const URL working_directory;
        Interface(const size_t argc, const char** argv);
        virtual ~Interface();
        inline const string& name() const {
            if(!command->name.empty()) {
                return command->name;
            } else {
                return application_name;
            }
        };
        inline Document operation() const {
            return selected->operation();
        };
        bool help_triggered() const {
            return selected->help_triggered();
        };
        bool version_triggered() const {
            return selected->version_triggered();
        };
        ostream& print_help(ostream& o) const;
        ostream& print_version(ostream& o) const;

    protected:
        Document configuration;
        Action* command;
        Action* selected;
        Layout layout;
        list< Action* > action_by_index;
        unordered_map< string, Action* > action_by_name;
        virtual ostream& print_action_dictionary(ostream& o) const;
        virtual void apply_action_base();
        virtual void load_action_array();
        virtual void load_selected_action();
        virtual void load_sub_action(const Value& ontology);
};

#endif /* PHENIQS_INTERFACE_H */
