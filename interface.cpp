
#include "interface.h"
#include "version.h"
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
    for(int i = 0; i < argc; i++) {
        value.append(argv[i]);
        value.append(" ");
    }
    return value;
};

void to_string(const ProgramAction& value, string& result) {
    switch(value) {
        case ProgramAction::DEMULTIPLEX:    result.assign("demux");      break;
        case ProgramAction::QUALITY:        result.assign("quality");    break;
        default:                            result.assign("unknown");    break;
    }
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

void to_string(const ParameterType& value, string& result) {
    switch(value) {
        case ParameterType::BOOLEAN:    result.assign("boolean");   break;
        case ParameterType::INTEGER:    result.assign("integer");   break;
        case ParameterType::DECIMAL:    result.assign("decimal");   break;
        case ParameterType::STRING:     result.assign("string");    break;
        case ParameterType::URL:        result.assign("url");       break;
        default:                        result.assign("unknown");   break;
    }
};
bool from_string(const char* value, ParameterType& result) {
         if(value == NULL)              result = ParameterType::UNKNOWN;
    else if(!strcmp(value, "boolean"))  result = ParameterType::BOOLEAN;
    else if(!strcmp(value, "integer"))  result = ParameterType::INTEGER;
    else if(!strcmp(value, "decimal"))  result = ParameterType::DECIMAL;
    else if(!strcmp(value, "string"))   result = ParameterType::STRING;
    else if(!strcmp(value, "url"))      result = ParameterType::URL;
    else                                result = ParameterType::UNKNOWN;
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

Prototype::Prototype(const Value& node) :
    cardinality(0),
    plural(false),
    mandatory(false),
    positional(false) {

    if(node.IsObject()) {
        decode_value_by_key< ParameterType >("type", type, node);
        decode_value_by_key< string >("name", name, node);
        decode_value_by_key< string >("help", help, node);
        decode_value_by_key< string >("meta", meta, node);
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
                default:
                    break;
            }
        }

        Value::ConstMemberIterator element;
        element = node.FindMember("cardinality");
        if(element != node.MemberEnd()) {
            if(element->value.IsString()) {
                if(element->value.GetStringLength() == 1 && *(element->value.GetString()) == '*') {
                    plural = true;
                } else { throw ConfigurationError("incorrect prototype syntax"); }
            } else if(element->value.IsUint64()) {
                plural = true;
                cardinality = element->value.GetUint64();
            } else { throw ConfigurationError("incorrect prototype syntax"); }
        }

        decode_value_by_key< bool >("mandatory", mandatory, node);

        element = node.FindMember("handle");
        if(element != node.MemberEnd()) {
            if(element->value.IsArray()) {
                for(auto& v : element->value.GetArray()) {
                    if(v.IsString()) {
                        SizeType length = v.GetStringLength();
                        const Value::Ch* value = v.GetString();
                        if(length == 2 && value[0] == '-' && value[1] != '-') {
                            // short
                            handles.emplace_back(value, 1, length - 1);
                        } else if(length > 3 && value[0] == '-' && value[1] == '-' && value[2] != '-') {
                            // long
                            handles.emplace_back(value, 2, length - 2);
                        } else { throw ConfigurationError("incorrect option handle syntax " + string(value, length)); }
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
             } else { throw ConfigurationError("incorrect prototype syntax"); }
        }
        positional = handles.empty();

        element = node.FindMember("choice");
        if(element != node.MemberEnd()) {
            if(element->value.IsArray()) {
                for(auto& v : element->value.GetArray()) {
                    if(v.IsString()) {
                        choices.emplace_back(v.GetString(), v.GetStringLength());
                    } else { throw ConfigurationError("incorrect prototype syntax"); }
                }
            } else { throw ConfigurationError("incorrect prototype syntax"); }
        }
    } else { throw ConfigurationError("incorrect prototype syntax"); }
};
Prototype::Prototype() :
    cardinality(0),
    plural(false),
    mandatory(false),
    positional(false) {
};
ostream& Prototype::print_help(ostream& o, const int& max_option_handle, const Layout& layout) const {
    o << setw(layout.option_indent) << ' ';
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
    o << setw(layout.complement_handle(max_option_handle, handle_length())) << ' ';
    o << help << endl;
    return o;
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
        length++;
        length += meta.length();
    }
    return length;
};
bool Prototype::is_choice() const {
    return choices.size() > 0;
};

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
            case ParameterType::URL: {
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
            case ParameterType::URL: {
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
            case ParameterType::URL: {
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
            case ParameterType::URL: {
                delete url_value;
                break;
            };
            default: {
                break;
            };
        };
    }
};
void Argument::set_value(const bool value) {
    assigned = true;
    *boolean_value = value;
};
void Argument::set_value(const int64_t value) {
    assigned = true;
    if(!prototype.plural) {
        *integer_value = value;
    } else {
        integer_array_value->push_back(value);
    }
};
void Argument::set_value(const double value) {
    assigned = true;
    if(!prototype.plural) {
        *decimal_value = value;
    } else {
        decimal_array_value->push_back(value);
    }
};
void Argument::set_value(const char* value) {
    assigned = true;
    if(!prototype.plural) {
        string_value->assign(value);
    } else {
        string_array_value->emplace_back(value);
    }
};
void Argument::set_value(const URL value) {
    assigned = true;
    if(!prototype.plural) {
        *url_value = value;
    } else {
        url_array_value->emplace_back(value);
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
uint64_t Argument::cardinality() const {
    uint64_t value(0);
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

Action::Action(const Value& node, bool root) :
    root(root),
    max_option_handle(0),
    default_value_node(NULL) {
    if(node.IsObject()) {
        decode_value_by_key< string >("name", name, node);
        decode_value_by_key< string >("description", description, node);
        decode_value_by_key< string >("epilog", epilog, node);

        Value::ConstMemberIterator element;
        element = node.FindMember("option");
        if(element != node.MemberEnd()) {
            if(element->value.IsArray()) {
                for(auto& o : element->value.GetArray()) {
                    Prototype* prototype = new Prototype(o);
                    if(!prototype->name.empty()) {
                        if(option_name_lookup.find(prototype->name) == option_name_lookup.end()) {
                            option_name_lookup[prototype->name] = prototype;
                            for(const auto& handle : prototype->handles) {
                                if(option_handle_lookup.find(handle) == option_handle_lookup.end()) {
                                    option_handle_lookup[handle] = prototype;
                                } else { throw ConfigurationError("redefining handle " + handle + " in option " + prototype->name); }
                            }
                            if(prototype->positional) {
                                positional_order.push_back(prototype);
                            } else {
                                optional_order.push_back(prototype);
                            }
                            max_option_handle = max(max_option_handle, prototype->handle_length());
                        } else { throw ConfigurationError("redefining " + prototype->name); }
                    } else { throw ConfigurationError("missing a name"); }
                }
            } else { throw ConfigurationError("incorrect syntax"); }
        }
        element = node.FindMember("default");
        if(element != node.MemberEnd()) {
            default_value_node = &element->value;
        }
    } else { throw ConfigurationError("incorrect action syntax"); }
};
Action::~Action() {
    for(auto prototype : optional_order) {
        delete prototype;
    }
    for(auto prototype : positional_order) {
        delete prototype;
    }
};
const Value* Action::default_value() {
    return default_value_node;
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
            case ParameterType::URL:
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
    for(const auto prototype : positional_order) {
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
                block.append(prototype->meta);
                break;
            case ParameterType::STRING: {
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
    if(!positional_order.empty()) {
        o << endl << endl << "Positional:" << endl;
        for(const auto option : positional_order) {
            option->print_help(o, max_option_handle, layout);
        }
    }
    if(!optional_order.empty()) {
        o << endl << endl << "Optional:" << endl;
        for(const auto option : optional_order) {
            option->print_help(o, max_option_handle, layout);
        }
    }
    return o;
};
bool Action::is_root() {
    return root;
};

CommandLine::CommandLine(const int argc, const char** argv) :
    argc(argc),
    argv(argv),
    application_name(argv[0]),
    application_version(PHENIQS_VERSION),
    full_command(assemble_full_command(argc, argv)),
    working_directory(get_cwd(), true),
    command(NULL),
    selected(NULL) {

    const string content(reinterpret_cast<const char *>(configuration_json), configuration_json_len);
    if(!_configuration.Parse(content.c_str()).HasParseError()) {
        if(_configuration.IsObject()) {
            encode_key_value("version", application_version, _configuration, _configuration);
        } else { throw ConfigurationError("interface configuration must be a dictionary"); }
    } else {
        throw ConfigurationError (
            string(GetParseError_En(_configuration.GetParseError())) + " at position " + to_string(_configuration.GetErrorOffset())
        );
    }

    command = new Action(_configuration, true);
    Value::ConstMemberIterator element = _configuration.FindMember("action");
    if(element != _configuration.MemberEnd()) {
        if(element->value.IsArray()) {
            for(auto& v : element->value.GetArray()) {
                Action* action = new Action(v, false);
                if(!action->name.empty()) {
                    if(action_name_lookup.find(action->name) == action_name_lookup.end()) {
                        action_name_lookup[action->name] = action;
                        action_order.push_back(action);
                        layout.max_action_name = max(layout.max_action_name, action->name_length());
                    } else { throw ConfigurationError("redefining " + action->name); }
                } else { throw ConfigurationError("action missing a name"); }
            }
        } else { throw ConfigurationError("incorrect action syntax"); }
    }
    decode();
    validate();
    load_instruction();
};
CommandLine::~CommandLine() {
    for(auto action : action_order) {
        delete action;
    }
    for(auto argument : argument_order) {
        delete argument;
    }
};
Document* CommandLine::load_default_instruction() {
    Document* document = new Document(kObjectType);
    if(default_value() != NULL) {
        document->CopyFrom(*default_value(), document->GetAllocator());
    }
    encode_key_value("working directory", working_directory, *document, *document);
    encode_key_value("base input url", working_directory, *document, *document);
    encode_key_value("base output url", working_directory, *document, *document);
    encode_key_value("action", selected->name, *document, *document);
    return document;
};
Document* CommandLine::load_instruction_from_configuration_file() {
    Document* document = NULL;
    URL* url(get_url("configuration url"));
    if(url != NULL) {
        if(url->is_readable()) {
            document = new Document();
            ifstream file(url->path());
            const string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
            file.close();
            if(document->Parse(content.c_str()).HasParseError()) {
                string message(GetParseError_En(document->GetParseError()));
                message += " at position ";
                message += to_string(document->GetErrorOffset());
                throw ConfigurationError(message);
            }
        } else { throw IOError("could not read configuration from " + string(*url)); }
    }
    return document;
};
Document* CommandLine::load_instruction_from_command_line(const Value& base) {
    Document* document = NULL;
    document = new Document(kObjectType);
    for(auto& record : base.GetObject()) {
        string key(record.name.GetString(), record.name.GetStringLength());
        Argument* argument = get(key);
        if(argument != NULL) {
            if(argument->prototype.plural) {
                switch(argument->prototype.type) {
                    case ParameterType::BOOLEAN: {
                        /* impossible */
                        break;
                    };
                    case ParameterType::INTEGER: {
                        list< int64_t >* value = get_integer_plural(key);
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::DECIMAL: {
                        list< double >* value = get_decimal_plural(key);
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::STRING: {
                        list< string >* value = get_string_plural(key);
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::URL: {
                        list< URL >* value = get_url_plural(key);
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    default:
                        break;
                }
            } else {
                switch(argument->prototype.type) {
                    case ParameterType::BOOLEAN: {
                        bool* value(get_boolean(key));
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::INTEGER: {
                        int64_t* value(get_integer(key));
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::DECIMAL: {
                        double* value(get_decimal(key));
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::STRING: {
                        string* value(get_string(key));
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    case ParameterType::URL: {
                        URL* value(get_url(key));
                        if(value != NULL) { encode_key_value(key, *value, *document, *document); }
                        break;
                    };
                    default:
                        break;
                }
            }
        }
    }
    return document;
};
void CommandLine::load_instruction() {
    Document* default_instruction = load_default_instruction();
    Document* instruction_from_file = load_instruction_from_configuration_file();
    Document* instruction_from_command_line = load_instruction_from_command_line(*default_instruction);

    if(instruction_from_file != NULL) {
        Document* base_configuration = new Document();
        merge_json_value(*default_instruction, *instruction_from_file, *base_configuration, *base_configuration);
        delete instruction_from_file;

        merge_json_value(*base_configuration, *instruction_from_command_line, _instruction, _instruction);
        delete base_configuration;

    } else {
        merge_json_value(*default_instruction, *instruction_from_command_line, _instruction, _instruction);

    }
    delete default_instruction;
    delete instruction_from_command_line;
};
const Document& CommandLine::instruction() const {
    return _instruction;
};
const string& CommandLine::name() const {
    if(!command->name.empty()) {
        return command->name;
    } else {
        return application_name;
    }
};
ProgramAction CommandLine::get_selected_action() const {
    ProgramAction value;
    from_string(selected->name, value);
    return value;
};
bool CommandLine::has_argument(const string& name) {
    return argument_name_lookup.count(name) > 0;
};
const Value* CommandLine::default_value() const {
    return selected->default_value();
};
Argument* CommandLine::get(const string& name) const {
    auto record = argument_name_lookup.find(name);
    if(record != argument_name_lookup.end()) {
        return record->second;
    } else {
        return NULL;
    }
};
bool* CommandLine::get_boolean(const string& name) const {
    bool* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_boolean();
    }
    return value;
};
int64_t* CommandLine::get_integer(const string& name) const {
    int64_t* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_integer();
    }
    return value;
};
double* CommandLine::get_decimal(const string& name) const {
    double* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_decimal();
    }
    return value;
};
string* CommandLine::get_string(const string& name) const {
    string* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_string();
    }
    return value;
};
URL* CommandLine::get_url(const string& name) const {
    URL* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_url();
    }
    return value;
};
list< int64_t >* CommandLine::get_integer_plural(const string& name) const {
    list< int64_t >* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_integer_array();
    }
    return value;
};
list< double >* CommandLine::get_decimal_plural(const string& name) const {
    list< double >* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_decimal_array();
    }
    return value;
};
list< string >* CommandLine::get_string_plural(const string& name) const {
    list< string >* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_string_array();
    }
    return value;
};
list< URL >* CommandLine::get_url_plural(const string& name) const {
    list< URL >* value(NULL);
    Argument* argument = get(name);
    if(argument) {
        value = argument->get_url_array();
    }
    return value;
};
bool CommandLine::help_triggered() const {
    if(argc < 2) {
        return true;
    } else {
        bool* value(get_boolean("help only"));
        if(value != NULL && value) {
            return true;
        } else {
            return false;
        }
    }
};
ostream& CommandLine::print_help(ostream& o) const {
    o << endl;
    print_version_element(o, layout);
    if(!selected->root) {
        command->print_description_element(o, layout);
    }
    selected->print_description_element(o, layout);
    selected->print_help(o, name(), layout);
    if(selected->root) {
        print_action_element(o, layout);
    }
    selected->print_epilog_element(o, layout);
    o << endl;
    return o;
};
bool CommandLine::version_triggered() const {
    bool* value(get_boolean("version only"));
    if(value != NULL && value) {
        return true;
    } else {
        return false;
    }
};
ostream& CommandLine::print_version(ostream& o) const {
    print_version_element(o, layout);
    return o;
};
Argument* CommandLine::parse_argument(const Prototype* prototype, size_t& index) {
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
                argument->set_value(argv[index]);
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
Argument* CommandLine::decode_optional(Action* action, size_t& index, const string& handle, bool& positional, bool composite) {
    Argument* argument(NULL);
    auto record = action->option_handle_lookup.find(handle);
    if(record != action->option_handle_lookup.end()) {
        Prototype* prototype = record->second;
        // if the option is not a boolean flag
        if(prototype->type != ParameterType::BOOLEAN) {
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
Argument* CommandLine::decode_positional(Action* action, size_t& index, const size_t& position) {
    Argument* argument = NULL;
    if(position < action->positional_order.size()) {
        Prototype* prototype = action->positional_order[position];
        argument = parse_argument(prototype, index);
    } else {
        throw CommandLineError("too many positional arguments");
    }
    return argument;
};
void CommandLine::decode() {
    size_t index(1);
    size_t position(0);
    bool positional(false);
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
                const char* key(argv[index]);
                size_t length(strlen(key));
                if(key[0] == '-' && length > 1) {
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
void CommandLine::validate() {
    if(!help_triggered()) {
        for(const auto& record : selected->option_name_lookup) {
            Argument* argument = NULL;
            const string& name = record.first;
            const Prototype* prototype = record.second;

            auto argument_record = argument_name_lookup.find(name);
            if(argument_record != argument_name_lookup.end()) {
                argument = argument_record->second;
            }
            if(prototype->mandatory && (argument == NULL || !argument->assigned)) {
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
void CommandLine::load_action(size_t& position) {
    if(argc > 1) {
        for(size_t i = 1; i < argc; i++) {
            if(argv[i][0] != '-') {
                string name(argv[i]);
                position = i + 1;
                auto record = action_name_lookup.find(name);
                if(record != action_name_lookup.end()) {
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
Argument* CommandLine::get_argument(const Prototype* prototype) {
    Argument* argument = NULL;
    auto record = argument_name_lookup.find(prototype->name);
    if(record != argument_name_lookup.end()) {
        argument = record->second;
    } else {
        argument = new Argument(prototype);
        argument_name_lookup[prototype->name] = argument;
        argument_order.push_back(argument);
    }
    return argument;
};
ostream& CommandLine::print_version_element(ostream& o, const Layout& layout) const {
    o << name() << " version " << application_version << endl;

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
ostream& CommandLine::print_action_element(ostream& o, const Layout& layout) const {
    o << endl << "available action" << endl;
    for(const auto action : action_order) {
        o << setw(layout.option_indent) << ' ';
        o << action->name;
        o << setw(layout.complement_action(action->name_length())) << ' ';
        o << action->description;
        o << endl;
    }
    return o;
};
