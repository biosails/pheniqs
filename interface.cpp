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
#include "configuration.h"

static inline string get_cwd() {
    char* buffer(NULL);
    char* temp(NULL);
    string directory;
    size_t size(128);

    if((buffer = static_cast< char* >(malloc(size))) == NULL) {
        throw OutOfMemoryError();
    }
    while(getcwd(buffer, size) == NULL) {
        switch(errno) {
            case EACCES: {
                free(buffer);
                throw IOError("insufficient permission to probe working directory");
                break;
            };
            case ERANGE: {
                size *= 2;
                if((temp = static_cast< char* >(realloc(buffer, size))) == NULL) {
                    free(buffer);
                    throw OutOfMemoryError();
                } else {
                    buffer = temp;
                }
                break;
            };
            default: {
                throw InternalError("error " + to_string(errno) + " when probing working directory");
                break;
            };
        }
    }
    directory.assign(buffer);
    free(buffer);
    return directory;
};
static inline string assemble_full_command(const int argc, const char** argv) {
    string value;
    value.append(argv[0]);
    for(int i = 1; i < argc; ++i) {
        value.append(" ");
        value.append(argv[i]);
    }
    return value;
};

void to_string(const ParameterType& value, string& result) {
    switch(value) {
        case ParameterType::BOOLEAN:    result.assign("boolean");   break;
        case ParameterType::INTEGER:    result.assign("integer");   break;
        case ParameterType::DECIMAL:    result.assign("decimal");   break;
        case ParameterType::STRING:     result.assign("string");    break;
        case ParameterType::URL:        result.assign("url");       break;
        case ParameterType::DIRECTORY:  result.assign("directory"); break;
        default:                        result.assign("unknown");   break;
    }
};
bool from_string(const char* value, ParameterType& result) {
         if(value == NULL)                  result = ParameterType::UNKNOWN;
    else if(!strcmp(value, "boolean"))      result = ParameterType::BOOLEAN;
    else if(!strcmp(value, "integer"))      result = ParameterType::INTEGER;
    else if(!strcmp(value, "decimal"))      result = ParameterType::DECIMAL;
    else if(!strcmp(value, "string"))       result = ParameterType::STRING;
    else if(!strcmp(value, "url"))          result = ParameterType::URL;
    else if(!strcmp(value, "directory"))    result = ParameterType::DIRECTORY;
    else                                    result = ParameterType::UNKNOWN;
    return (result == ParameterType::UNKNOWN ? false : true);
};
bool from_string(const string& value, ParameterType& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const ParameterType& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const ParameterType& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< ParameterType >(const Value::Ch* key, ParameterType& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template<> ParameterType decode_value_by_key< ParameterType >(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsString()) {
                    ParameterType value;
                    from_string(reference->value.GetString(), value);
                    return value;
                } else { throw ConfigurationError(string(key) + " element must be a string"); }
            } else { throw ConfigurationError(string(key) + " is null"); }
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};

/*  Prototype */

Prototype::Prototype(const Value& ontology) :
    type(decode_value_by_key< ParameterType >("type", ontology)),
    name(decode_value_by_key< string >("name", ontology)),
    cardinality(0),
    plural(false),
    mandatory(false),
    positional(false) {

    decode_value_by_key< string >("help", help, ontology);
    decode_value_by_key< bool >("mandatory", mandatory, ontology);
    if(!decode_value_by_key< string >("meta", meta, ontology)) {
        if(meta.empty()) {
            switch(type) {
                case ParameterType::BOOLEAN:
                    meta.clear();
                    break;
                case ParameterType::INTEGER:
                    meta.assign("INT");
                    break;
                case ParameterType::DECIMAL:
                    meta.assign("FLOAT");
                    break;
                case ParameterType::STRING:
                    meta.assign("STRING");
                    break;
                case ParameterType::URL:
                    meta.assign("URL");
                    break;
                case ParameterType::DIRECTORY:
                    meta.assign("DIR");
                    break;
                default:
                    break;
            }
        }
    }

    Value::ConstMemberIterator reference = ontology.FindMember("cardinality");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsString()) {
            if(reference->value.GetStringLength() == 1 && *(reference->value.GetString()) == '*') {
                plural = true;
            } else { throw ConfigurationError("option cardinality must be an intger or *"); }
        } else if(reference->value.IsUint()) {
            plural = true;
            cardinality = reference->value.GetUint();
        } else { throw ConfigurationError("option cardinality must be an intger or *"); }
    }

    reference = ontology.FindMember("handle");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsArray()) {
            for(auto& v : reference->value.GetArray()) {
                if(v.IsString()) {
                    SizeType length = v.GetStringLength();
                    const Value::Ch* value = v.GetString();
                    if(length == 2 && value[0] == '-' && value[1] != '-') {
                        /* short bsd option syntax */
                        handles.emplace_back(value, 1, length - 1);
                    } else if(length > 3 && value[0] == '-' && value[1] == '-' && value[2] != '-') {
                        /* long bsd option syntax */
                        handles.emplace_back(value, 2, length - 2);
                    } else { throw ConfigurationError("incorrect option handle syntax " + string(value, length)); }
                } else { throw ConfigurationError("option handle must be a string"); }
            }
        } else { throw ConfigurationError("option handle element must be an array"); }
    }
    positional = handles.empty();

    reference = ontology.FindMember("choice");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsArray()) {
            for(auto& v : reference->value.GetArray()) {
                if(v.IsString()) {
                    choices.emplace_back(v.GetString(), v.GetStringLength());
                } else { throw ConfigurationError("option choice array element must be a string"); }
            }
        } else { throw ConfigurationError("option choice element must be an array"); }
    }
};
int Prototype::handle_length() const {
    int length(0);
    if(!positional) {
        for(const auto& handle : handles) {
            length += handle.length();
            length += (handle.length() > 1 ? 4 : 3);
        }
        length--;
    }
    if(meta.length() > 0) {
        ++length;
        length += meta.length();
    }
    return length;
};
ostream& Prototype::print_help(ostream& o, const int& max_option_handle, const Layout& layout) const {
    o << setw(layout.option_indent) << ' ';
    if(!positional) {
        for(size_t i = 0; i < handles.size(); ++i) {
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
    o << setw(layout.complement_handle(max_option_handle, handle_length())) << ' ';
    o << help << endl;
    return o;
};

/*  Argument */

Argument::Argument(const Prototype* prototype) :
    assigned(false),
    prototype(*prototype) {
    if(prototype->plural) {
        switch(prototype->type) {
            case ParameterType::INTEGER: {
                integer_array_value = new list< int64_t >();
                break;
            };
            case ParameterType::DECIMAL: {
                decimal_array_value = new list< double >();
                break;
            };
            case ParameterType::STRING: {
                string_array_value = new list< string >();
                break;
            };
            case ParameterType::URL:
            case ParameterType::DIRECTORY: {
                url_array_value = new list< URL >();
                break;
            };
            default: {
                break;
            };
        };
    } else {
        switch(prototype->type) {
            case ParameterType::BOOLEAN: {
                boolean_value = new bool(false);
                break;
            };
            case ParameterType::INTEGER: {
                integer_value = new int64_t(0);
                break;
            };
            case ParameterType::DECIMAL: {
                decimal_value = new double(0);
                break;
            };
            case ParameterType::STRING: {
                string_value = new string();
                break;
            };
            case ParameterType::URL:
            case ParameterType::DIRECTORY: {
                url_value = new URL();
                break;
            };
            default: {
                break;
            };
        };
    }
};
Argument::~Argument() {
    if(prototype.plural) {
        switch(prototype.type) {
            case ParameterType::INTEGER: {
                delete integer_array_value;
                break;
            };
            case ParameterType::DECIMAL: {
                delete decimal_array_value;
                break;
            };
            case ParameterType::STRING: {
                delete string_array_value;
                break;
            };
            case ParameterType::URL:
            case ParameterType::DIRECTORY: {
                delete url_array_value;
                break;
            };
            default: {
                break;
            };
        };
    } else {
        switch(prototype.type) {
            case ParameterType::BOOLEAN: {
                delete boolean_value;
                break;
            };
            case ParameterType::INTEGER: {
                delete integer_value;
                break;
            };
            case ParameterType::DECIMAL: {
                delete decimal_value;
                break;
            };
            case ParameterType::STRING: {
                delete string_value;
                break;
            };
            case ParameterType::URL:
            case ParameterType::DIRECTORY: {
                delete url_value;
                break;
            };
            default: {
                break;
            };
        };
    }
};
bool* Argument::get_boolean() const {
    bool* value(NULL);
    if(prototype.type == ParameterType::BOOLEAN) {
        value = boolean_value;
    }
    return value;
};
int64_t* Argument::get_integer() const {
    int64_t* value(NULL);
    if(prototype.type == ParameterType::INTEGER) {
        value = integer_value;
    }
    return value;
};
double* Argument::get_decimal() const {
    double* value(NULL);
    if(prototype.type == ParameterType::DECIMAL) {
        value = decimal_value;
    }
    return value;
};
string* Argument::get_string() const {
    string* value(NULL);
    if(prototype.type == ParameterType::STRING) {
        value = string_value;
    }
    return value;
};
URL* Argument::get_url() const {
    URL* value(NULL);
    if(prototype.type == ParameterType::URL) {
        value = url_value;
    }
    return value;
};
URL* Argument::get_directory() const {
    URL* value(NULL);
    if(prototype.type == ParameterType::DIRECTORY) {
        value = url_value;
    }
    return value;
};
list< int64_t >* Argument::get_integer_array() const {
    list< int64_t >* value(NULL);
    if(prototype.type == ParameterType::INTEGER && prototype.plural) {
        value = integer_array_value;
    }
    return value;
};
list< double >* Argument::get_decimal_array() const {
    list< double >* value(NULL);
    if(prototype.type == ParameterType::DECIMAL && prototype.plural) {
        value = decimal_array_value;
    }
    return value;
};
list< string >* Argument::get_string_array() const {
    list< string >* value(NULL);
    if(prototype.type == ParameterType::STRING && prototype.plural) {
        value = string_array_value;
    }
    return value;
};
list< URL >* Argument::get_url_array() const {
    list< URL >* value(NULL);
    if(prototype.type == ParameterType::URL && prototype.plural) {
        value = url_array_value;
    }
    return value;
};
list< URL >* Argument::get_directory_array() const {
    list< URL >* value(NULL);
    if(prototype.type == ParameterType::DIRECTORY && prototype.plural) {
        value = url_array_value;
    }
    return value;
};
uint32_t Argument::cardinality() const {
    uint32_t value(0);
    if(!prototype.plural) {
        value = assigned ? 1 : 0;
    } else {
        switch(prototype.type) {
            case ParameterType::BOOLEAN:
                // does not make sense
                break;
            case ParameterType::INTEGER:
                value = integer_array_value->size();
                break;
            case ParameterType::DECIMAL:
                value = decimal_array_value->size();
                break;
            case ParameterType::STRING:
                value = string_array_value->size();
                break;
            case ParameterType::URL:
            case ParameterType::DIRECTORY:
                value = url_array_value->size();
                break;
            default:
                break;
        }
    }
    return value;
};
bool Argument::satisfied() const {
    bool value(true);
    if(prototype.plural) {
        if(prototype.cardinality == 0 || prototype.cardinality > cardinality()) {
            value = false;
        }
    } else {
        if(!assigned) {
            value = false;
        }
    }
    return value;
};

/*  Action */

Action::Action(const Value& ontology, bool root) :
    ontology(ontology),
    name(decode_value_by_key< string >("name", ontology)),
    description(decode_value_by_key< string >("description", ontology)),
    epilog(decode_value_by_key< string >("epilog", ontology)),
    root(root),
    max_option_handle(0) {

    Value::ConstMemberIterator reference = ontology.FindMember("option");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.IsArray()) {
            for(auto& element : reference->value.GetArray()) {
                if(element.IsObject()) {
                    string key;
                    if(decode_value_by_key< string >("name", key, element)) {
                        Prototype* prototype = new Prototype(element);
                        if(option_by_name.find(key) == option_by_name.end()) {
                            option_by_name.emplace(make_pair(key, prototype));
                            for(const auto& handle : prototype->handles) {
                                if(option_by_handle.find(handle) == option_by_handle.end()) {
                                    option_by_handle[handle] = prototype;
                                } else { throw ConfigurationError("redefining handle " + handle + " in option " + prototype->name); }
                            }
                            if(prototype->positional) {
                                positional_by_index.push_back(prototype);
                            } else {
                                optional_by_index.push_back(prototype);
                            }
                            max_option_handle = max(max_option_handle, prototype->handle_length());
                        } else { throw ConfigurationError("redefining option " + prototype->name); }
                    }
                } else { throw ConfigurationError("option element must be a dictionary"); }
            }
        } else { throw ConfigurationError("incorrect syntax"); }
    }
};
Action::~Action() {
    for(auto prototype : optional_by_index) {
        delete prototype;
    }
    for(auto prototype : positional_by_index) {
        delete prototype;
    }
    for(auto argument : argument_by_index) {
        delete argument;
    }
};
Argument* Action::get_argument(const Prototype* prototype) {
    Argument* argument = NULL;
    auto record = argument_by_name.find(prototype->name);
    if(record != argument_by_name.end()) {
        argument = record->second;
    } else {
        argument = new Argument(prototype);
        argument_by_name[prototype->name] = argument;
        argument_by_index.push_back(argument);
    }
    return argument;
};
Argument* Action::parse_argument(const size_t argc, const char** argv, const Prototype* prototype, size_t& index) {
    Argument* argument(NULL);
    try {
        switch(prototype->type) {
            case ParameterType::BOOLEAN: {
                argument = get_argument(prototype);
                argument->set_value(true);
                break;
            };
            case ParameterType::INTEGER: {
                string raw(argv[index]);
                argument = get_argument(prototype);
                argument->set_value(static_cast < int64_t>(stoll(raw)));
                break;
            };
            case ParameterType::DECIMAL: {
                string raw(argv[index]);
                argument = get_argument(prototype);
                argument->set_value(stod(raw));
                break;
            };
            case ParameterType::STRING: {
                argument = get_argument(prototype);
                argument->set_value(argv[index]);
                break;
            };
            case ParameterType::URL: {
                argument = get_argument(prototype);
                string buffer(argv[index]);
                URL url(expand_shell(buffer));
                argument->set_value(url);
                break;
            };
            case ParameterType::DIRECTORY: {
                argument = get_argument(prototype);
                string buffer(argv[index]);
                URL url(expand_shell(buffer), true);
                argument->set_value(url);
                break;
            };
            default:
                break;
        }
    } catch(invalid_argument& e) {
        throw CommandLineError("invalid value for argument " + prototype->name);
    } catch(out_of_range& e) {
        throw CommandLineError("out of range value for argument " + prototype->name);
    }
    return argument;
};
Argument* Action::load_optional_argument(const size_t argc, const char** argv, size_t& index, const string& handle, bool& positional, bool composite) {
    Argument* argument(NULL);
    auto record = option_by_handle.find(handle);
    if(record != option_by_handle.end()) {
        Prototype* prototype = record->second;
        // if the option is not a boolean flag
        if(prototype->type != ParameterType::BOOLEAN) {
            if(composite) {
                // if the handle was grouped and not the last in the group it must be a boolean option
                throw CommandLineError("argument " + prototype->name + " requires a value");
            } else  if(index + 1 < argc) {
                // if there is no next token or the next token is a noop than the option is missing a value
                ++index;
                if(!strcmp(argv[index], "--")) {
                    positional = true;
                    ++index;
                    throw CommandLineError("argument " + prototype->name + " requires a value");
                }
            } else { throw CommandLineError("argument " + prototype->name + " requires a value"); }
        }
        argument = parse_argument(argc, argv, prototype, index);
    }
    return argument;
};
Argument* Action::load_positional_argument(const size_t argc, const char** argv, size_t& index, const size_t& position) {
    Argument* argument = NULL;
    if(position < positional_by_index.size()) {
        Prototype* prototype = positional_by_index[position];
        argument = parse_argument(argc, argv, prototype, index);
    } else {
        throw CommandLineError("too many positional arguments");
    }
    return argument;
};
void Action::load(const size_t argc, const char** argv) {
    parse_command_line(argc, argv);
    validate_command_line();
};
void Action::parse_command_line(const size_t argc, const char** argv) {
    /*  if this is called on a sub command argc and argv are already stripped of the main command
        so argv[0] is the sub command name. This makes this gives it the same behaviour as calling
        it on the main command without a sub command.
    */
    if(argc > 1) {
        size_t index(1);
        size_t position(0);
        bool positional(false);

        while(index < argc) {
            if(!strcmp(argv[index], "--")) {
                /* break signal encountered, this means no more optional arguments */
                if(positional) {
                    /* if already parsing positional arguments skip to the next positional */
                    ++position;
                } else {
                    /* if currently parsing optional arguments switch to positional parsing */
                    positional = true;
                }
            } else {
                if(!positional) {
                    const char* key(argv[index]);
                    size_t length(strlen(key));
                    if(key[0] == '-' && length > 1) {
                        ++key;
                        if(*key == '-') {
                            if(length > 3) {
                                /* long bsd syntax with -- prefix */
                                ++key;
                                string handle(key);
                                if(!load_optional_argument(argc, argv, index, handle, positional)) {
                                    throw CommandLineError("unknown argument --" + handle);
                                }
                            } else {
                                throw CommandLineError("unknown argument " + string(argv[index]));
                            }
                        } else {
                            /* short bsd syntax, potentially clustered */
                            while(*key != '\0') {
                                string handle(1, *key);
                                if(!load_optional_argument(argc, argv, index, handle, positional, *(key + 1) != '\0')) {
                                    throw CommandLineError("unknown argument -" + handle);
                                }
                                ++key;
                            }
                        }
                    } else {
                        /* first positional argument encountered, switch to positional parsing */
                        positional = true;
                    }
                }
                if(positional) {
                    if(position < positional_by_index.size()) {
                        Argument* argument = load_positional_argument(argc, argv, index, position);
                        if(argument->satisfied()) {
                            ++position;
                        }
                    }
                }
            }
            ++index;
        }
    } else { set_help_triggered(); }
};
void Action::set_help_triggered() {
    auto record = option_by_handle.find("help");
    if(record != option_by_handle.end()) {
        Prototype* prototype = record->second;
        Argument* argument = get_argument(prototype);
        argument->set_value(true);
    }
};
void Action::validate_command_line() {
    if(!help_triggered() && !version_triggered()) {
        for(const auto& record : option_by_name) {
            Argument* argument = NULL;
            const string& name = record.first;
            const Prototype* prototype = record.second;

            auto argument_record = argument_by_name.find(name);
            if(argument_record != argument_by_name.end()) {
                argument = argument_record->second;
            }
            if(prototype->mandatory && (argument == NULL || argument->empty())) {
                throw CommandLineError("missing mandatory argument " + prototype->name);
            } else if(prototype->plural && prototype->cardinality > 0 && (argument == NULL || prototype->cardinality != argument->cardinality())) {
                throw CommandLineError("argument " + prototype->name + " requires exactly " + to_string(prototype->cardinality) + " values");
            } else if(argument != NULL && prototype->is_choice()) {
                switch(prototype->type) {
                    case ParameterType::BOOLEAN: {
                        break;
                    };
                    case ParameterType::INTEGER: {
                        break;
                    };
                    case ParameterType::DECIMAL: {
                        break;
                    };
                    case ParameterType::STRING: {
                        bool found(false);
                        string* value(argument->get_string());
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
                    case ParameterType::URL: {
                        break;
                    };
                    default:
                        break;
                };
            }
        }
    }
};
ostream& Action::print_usage(ostream& o, const string& application_name, const Layout& layout) const {
    string buffer;
    string block;
    buffer.append("Usage : ");
    buffer.append(application_name);
    if(!root) {
        buffer.append(" ");
        buffer.append(name);
    }
    int indent(static_cast< int >(buffer.length()));
    for(const auto prototype : optional_by_index) {
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
            case ParameterType::BOOLEAN:
                block.append(prototype->handles[0]);
                break;
            case ParameterType::INTEGER:
                block.append(prototype->handles[0]);
                block.append(" ");
                block.append(prototype->meta);
                break;
            case ParameterType::DECIMAL:
                block.append(prototype->handles[0]);
                block.append(" ");
                block.append(prototype->meta);
                break;
            case ParameterType::STRING:
                block.append(prototype->handles[0]);
                block.append(" ");
                if(prototype->is_choice()) {
                    for(size_t i = 0; i < prototype->choices.size(); ++i) {
                        if(i > 0) {
                            block.append("|");
                        }
                        block.append(prototype->choices[i]);
                    }
                } else {
                    block.append(prototype->meta);
                }
                break;
            case ParameterType::URL:
            case ParameterType::DIRECTORY:
                block.append(prototype->handles[0]);
                block.append(" ");
                block.append(prototype->meta);
                break;
            default:
                break;
        }
        if(!prototype->mandatory) {
            block.append("]");
        }
        if(buffer.length() + block.length() > layout.max_line_width) {
            o << buffer << endl << setw(indent) << ' ';
            buffer.clear();
        }
        buffer.append(block);
        block.clear();
    }
    for(const auto prototype : positional_by_index) {
        block.append(" ");
        if(!prototype->mandatory) {
            block.append("[");
        }
        switch(prototype->type) {
            case ParameterType::BOOLEAN:
                break;
            case ParameterType::INTEGER:
            case ParameterType::DECIMAL:
            case ParameterType::URL:
            case ParameterType::DIRECTORY:
                block.append(prototype->meta);
                break;
            case ParameterType::STRING: {
                if(prototype->is_choice()) {
                    for(size_t i = 0; i < prototype->choices.size(); ++i) {
                        if(i > 0) { block.append("|"); }
                        block.append(prototype->choices[i]);
                    }
                } else {
                    block.append(prototype->meta);
                }
                break;
            };
            default:
                break;
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
        if(buffer.length() + block.length() > layout.max_line_width) {
            o << buffer << endl << setw(indent) << ' ';
            buffer.clear();
        }
        buffer.append(block);
        block.clear();
    }
    if(root) {
        buffer.append(" ACTION ...");
        if(buffer.length() + block.length() > layout.max_line_width) {
            o << buffer << endl << setw(indent) << ' ';
            buffer.clear();
        }
        buffer.append(block);
        block.clear();
    }
    o << buffer;
    return o;
};
ostream& Action::print_description_element(ostream& o, const Layout& layout) const {
    if(description.length() > 0) {
        o << description << endl << endl;
    }
    return o;
};
ostream& Action::print_epilog_element(ostream& o, const Layout& layout) const {
    if(epilog.length() > 0) {
        o << endl << epilog << endl;
    }
    return o;
};
ostream& Action::print_help(ostream& o, const string& application_name, const Layout& layout) const {
    print_usage(o, application_name, layout);
    if(!positional_by_index.empty()) {
        o << endl << endl << "Positional:" << endl;
        for(const auto option : positional_by_index) {
            option->print_help(o, max_option_handle, layout);
        }
    }
    if(!optional_by_index.empty()) {
        o << endl << endl << "Optional:" << endl;
        for(const auto option : optional_by_index) {
            option->print_help(o, max_option_handle, layout);
        }
    }
    return o;
};
Document Action::operation() {
    Document document(kObjectType);
    document.CopyFrom(ontology, document.GetAllocator());
    document.RemoveMember("option");

    Value interactive(kObjectType);
    for(auto& record : option_by_name) {
        const string& key = record.first;
        Argument* argument = get(key);
        if(argument != NULL) {
            if(argument->plural()) {
                switch(argument->type()) {
                    case ParameterType::BOOLEAN: {
                        /* impossible */
                        break;
                    };
                    case ParameterType::INTEGER: {
                        list< int64_t >* reference = get_integer_plural(key);
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::DECIMAL: {
                        list< double >* reference = get_decimal_plural(key);
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::STRING: {
                        list< string >* reference = get_string_plural(key);
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::URL:
                    case ParameterType::DIRECTORY: {
                        list< URL >* reference = get_url_plural(key);
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    default:
                        break;
                }
            } else {
                switch(argument->type()) {
                    case ParameterType::BOOLEAN: {
                        bool* reference(get_boolean(key));
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::INTEGER: {
                        int64_t* reference(get_integer(key));
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::DECIMAL: {
                        double* reference(get_decimal(key));
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::STRING: {
                        string* reference(get_string(key));
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    case ParameterType::URL:
                    case ParameterType::DIRECTORY: {
                        URL* reference(get_url(key));
                        if(reference != NULL) { encode_key_value(key, *reference, interactive, document); }
                        break;
                    };
                    default:
                        break;
                }
            }
        }
    }
    document.AddMember("interactive", interactive.Move(), document.GetAllocator());

    sort_json_value(document, document);
    return document;
};

/*  Interface */

Interface::Interface(const size_t argc, const char** argv) :
    argc(argc),
    argv(argv),
    application_name(argv[0]),
    application_version(PHENIQS_VERSION),
    full_command(assemble_full_command(argc, argv)),
    working_directory(get_cwd(), true),
    command(NULL),
    selected(NULL) {

    /* load the interface configuration */
    if(!configuration.Parse(configuration_json, configuration_json_len).HasParseError()) {
        if(configuration.IsObject()) {
            apply_action_base();
            load_action_array();
            load_selected_action();
        } else { throw ConfigurationError("interface configuration must be a dictionary"); }
    } else { throw ConfigurationError(string(GetParseError_En(configuration.GetParseError())) + " at position " + to_string(configuration.GetErrorOffset())); }
};
Interface::~Interface() {
    for(auto action : action_by_index) {
        delete action;
    }
};
void Interface::apply_action_base() {
    /* get the base code and decoder node from the configuration */
    Value default_configuration_action(kObjectType);
    Value::MemberIterator reference = configuration.FindMember("projection");
    if(reference != configuration.MemberEnd()) {
        Value& base = reference->value;
        if(!base.IsNull()) {
            reference = base.FindMember("action");
            if(reference != base.MemberEnd() && !reference->value.IsNull()) {
                default_configuration_action.CopyFrom(reference->value, configuration.GetAllocator());
            }
        }
    }

    reference = configuration.FindMember("default");
    if(reference != configuration.MemberEnd()) {
        if(!reference->value.IsNull()) {
            encode_key_value("working directory", working_directory, reference->value, configuration);
            encode_key_value("base input url", working_directory, reference->value, configuration);
            encode_key_value("base output url", working_directory, reference->value, configuration);
            encode_key_value("application version", application_version, reference->value, configuration);
            encode_key_value("application name", application_name, reference->value, configuration);
            encode_key_value("full command", full_command, reference->value, configuration);
        }
    }

    encode_key_value("working directory", working_directory, default_configuration_action, configuration);
    encode_key_value("base input url", working_directory, default_configuration_action, configuration);
    encode_key_value("base output url", working_directory, default_configuration_action, configuration);
    encode_key_value("application version", application_version, default_configuration_action, configuration);
    encode_key_value("application name", application_name, default_configuration_action, configuration);
    encode_key_value("full command", full_command, default_configuration_action, configuration);

    /* project action attributes from the document root */
    Value action_template(kObjectType);
    project_json_value(default_configuration_action, configuration, action_template, configuration);
    default_configuration_action.SetNull();

    /* no need for the action projection in the action configuration nodes */
    reference = action_template.FindMember("projection");
    if(reference != action_template.MemberEnd()) {
        if(!reference->value.IsNull()) {
            reference->value.RemoveMember("action");
        }
    }

    /* apply the action template on the root action */
    merge_json_value(action_template, configuration, configuration);

    /* apply the action template on the sub actions */
    reference = configuration.FindMember("action");
    if(reference != configuration.MemberEnd()) {
        if(reference->value.IsArray()) {
            for(auto& element : reference->value.GetArray()) {
                if(element.IsObject()) {
                    merge_json_value(action_template, element, configuration);
                } else { throw ConfigurationError("action element must be a dictionary"); }
            }
        } else { throw ConfigurationError("interface action element must be an array"); }
    }

    /* clean up template object */
    action_template.SetNull();
};
void Interface::load_action_array() {
    /* load the root action */
    command = new Action(configuration, true);

    /* load sub action */
    Value::MemberIterator reference = configuration.FindMember("action");
    if(reference != configuration.MemberEnd()) {
        if(reference->value.IsArray()) {
            for(auto& element : reference->value.GetArray()) {
                load_sub_action(element);
            }
        } else { throw ConfigurationError("interface action element must be an array"); }
    }
};
void Interface::load_sub_action(const Value& ontology) {
    string key;
    if(decode_value_by_key< string >("name", key, ontology)) {
        if(action_by_name.find(key) == action_by_name.end()) {
            Action* action = new Action(ontology);
            action_by_index.push_back(action);
            action_by_name.emplace(make_pair(key, action));
            layout.max_action_name = max(layout.max_action_name, action->name_length());
        } else { throw ConfigurationError("duplicate action " + key); }
    } else { throw InternalError("interface action must declare a name"); }
};
void Interface::load_selected_action() {
    /* start by assuming there is no sub command and the action is the main command */
    selected = command;
    size_t action_argc(argc);
    const char** action_argv(argv);

    /* identify the sub command  */
    if(!action_by_index.empty() && argc > 1) {
        /*  Look for the first argv element that is not a bsd option starting with  - */
        for(size_t i = 1; i < argc; ++i) {
            if(argv[i][0] != '-') {
                string key(argv[i]);
                auto record = action_by_name.find(key);
                if(record != action_by_name.end()) {
                    selected = record->second;
                    action_argc -= i;
                    action_argv += i;
                } else { throw CommandLineError("unknown action " +  key); }
                break;
            }
        }
    }
    selected->load(action_argc, action_argv);
};
ostream& Interface::print_help(ostream& o) const {
    o << endl;
    print_version(o);
    if(!selected->is_root()) {
        command->print_description_element(o, layout);
    }
    selected->print_description_element(o, layout);
    selected->print_help(o, name(), layout);
    if(selected->is_root()) {
        print_action_element(o);
    }
    selected->print_epilog_element(o, layout);
    o << endl;
    return o;
};
ostream& Interface::print_version(ostream& o) const {
    o << name() << " version " << application_version << endl;
    return o;
};
ostream& Interface::print_action_element(ostream& o) const {
    o << endl << "available action" << endl;
    for(const auto action : action_by_index) {
        o << setw(layout.option_indent) << ' ';
        o << action->name;
        o << setw(layout.complement_action(action->name_length())) << ' ';
        o << action->description;
        o << endl;
    }
    return o;
};
