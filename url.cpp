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

#include "url.h"

void to_string(const FormatType& value, string& result) {
    switch (value) {
        case FormatType::FASTQ: result.assign("fastq");  break;
        case FormatType::SAM:   result.assign("sam");    break;
        case FormatType::BAM:   result.assign("bam");    break;
        case FormatType::BAI:   result.assign("bai");    break;
        case FormatType::CRAM:  result.assign("cram");   break;
        case FormatType::CRAI:  result.assign("crai");   break;
        case FormatType::VCF:   result.assign("vcf");    break;
        case FormatType::BCF:   result.assign("bcf");    break;
        case FormatType::CSI:   result.assign("csi");    break;
        case FormatType::GZI:   result.assign("gzi");    break;
        case FormatType::TBI:   result.assign("tbi");    break;
        case FormatType::BED:   result.assign("bed");    break;
        case FormatType::JSON:  result.assign("json");   break;
        default:                                         break;
    }
};
bool from_string(const char* value, FormatType& result) {
         if(value == NULL)              result = FormatType::UNKNOWN;
    else if(!strcmp(value, "fastq"))    result = FormatType::FASTQ;
    else if(!strcmp(value, "sam"))      result = FormatType::SAM;
    else if(!strcmp(value, "bam"))      result = FormatType::BAM;
    else if(!strcmp(value, "bai"))      result = FormatType::BAI;
    else if(!strcmp(value, "cram"))     result = FormatType::CRAM;
    else if(!strcmp(value, "crai"))     result = FormatType::CRAI;
    else if(!strcmp(value, "vcf"))      result = FormatType::VCF;
    else if(!strcmp(value, "bcf"))      result = FormatType::BCF;
    else if(!strcmp(value, "csi"))      result = FormatType::CSI;
    else if(!strcmp(value, "gzi"))      result = FormatType::GZI;
    else if(!strcmp(value, "TBI"))      result = FormatType::TBI;
    else if(!strcmp(value, "bed"))      result = FormatType::BED;
    else if(!strcmp(value, "json"))     result = FormatType::JSON;
    else                                result = FormatType::UNKNOWN;
    return (result == FormatType::UNKNOWN ? false : true);
};
bool from_string(const string& value, FormatType& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const FormatType& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const FormatType& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< FormatType >(const Value::Ch* key, FormatType& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> FormatType decode_value_by_key(const Value::Ch* key, const Value& container) {
    FormatType value(FormatType::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};



void to_string(const IoDirection& value, string& result) {
    switch (value) {
        case IoDirection::IN:   result.assign("in");    break;
        case IoDirection::OUT:  result.assign("out");   break;
        default:                                        break;
    }
};
bool from_string(const char* value, IoDirection& result) {
         if(value == NULL)          result = IoDirection::UNKNOWN;
    else if(!strcmp(value, "in"))   result = IoDirection::IN;
    else if(!strcmp(value, "out"))  result = IoDirection::OUT;
    else                            result = IoDirection::UNKNOWN;
    return (result == IoDirection::UNKNOWN ? false : true);
};
bool from_string(const string& value, IoDirection& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const IoDirection& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const IoDirection& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< IoDirection >(const Value::Ch* key, IoDirection& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> IoDirection decode_value_by_key(const Value::Ch* key, const Value& container) {
    IoDirection value(IoDirection::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

URL::URL() :
    _type(FormatType::UNKNOWN) {
};
URL::URL(const URL& other) :
    _path(other._path),
    _basename(other._basename),
    _dirname(other._dirname),
    _extension(other._extension),
    _compression(other._compression),
    _type(other._type) {
};
URL::URL(const string& path) :
    _type(FormatType::UNKNOWN) {
    parse_file(path, IoDirection::UNKNOWN);
};
URL::URL(const string& path, const bool& is_directory) :
    _type(FormatType::UNKNOWN) {
    if(is_directory) {
        parse_directory(path);
    } else {
        parse_file(path, IoDirection::UNKNOWN);
    }
};
URL::URL(const string& path, const IoDirection& direction) :
    _type(FormatType::UNKNOWN) {
    parse_file(path, direction);
};
void URL::parse_directory(const string& path) {
    clear();
    if(!path.empty()) {
        _path.assign(path);

        // first resolve any environment variables in the path
        expand(_path);

        /* trim trailing file separator if present */
        if(_path.size() > 1 && _path.back() == PATH_SEPARATOR) {
            _path.erase(_path.size() - 1);
        }
        _dirname.assign(_path);
    }
};
void URL::parse_file(const string& path, const IoDirection& direction) {
    clear();
    if(!path.empty()) {
        _path.assign(path);

        // first resolve any environment variables in the path
        expand(_path);

        /*  standard stream handling
            first allow - as an alias for stdin and stdout according to IO direction
            otherwise check for non canonical paths to standard streams and replace with
            the canonical one
        */
        if(_path == STANDARD_STREAM_ALIAS) {
            switch(direction) {
                case IoDirection::IN: {
                    _path.assign(CANONICAL_STDIN_PATH);
                    break;
                };
                case IoDirection::OUT: {
                    _path.assign(CANONICAL_STDOUT_PATH);
                    break;
                };
                case IoDirection::UNKNOWN: {
                    throw ConfigurationError("can not interpret standard stream alias - without a declared IO direction");
                    break;
                };
            }
        } else {
            if(_path == "/dev/stdin" || _path == "/dev/fd/0" || _path == "/proc/self/fd/0") {
                if(direction == IoDirection::OUT) {
                    throw ConfigurationError("can not use " + _path + " for output");
                }
                _path.assign(CANONICAL_STDIN_PATH);

            } else if(_path == "/dev/stdout" || _path == "/dev/fd/1" || _path == "/proc/self/fd/1") {
                if(direction == IoDirection::IN) {
                    throw ConfigurationError("can not use " + _path + " for input");
                }
                _path.assign(CANONICAL_STDOUT_PATH);

            } else if(_path == "/dev/stderr" || _path == "/dev/fd/2" || _path == "/proc/self/fd/2") {
                if(direction == IoDirection::IN) {
                    throw ConfigurationError("can not use " + _path + " for input");
                }
                _path.assign(CANONICAL_STDERR_PATH);
            }
        }

        // split the path into dirname and basename
        auto position = _path.find_last_of(PATH_SEPARATOR);
        if(position != string::npos) {
            // the path separator is present
            if(position  + 1 < _path.size()) {
                _basename.assign(_path, position + 1, string::npos);
            }
            if(position > 0) {
                // the last path separator is not the first character
                _dirname.assign(_path, 0, position);
            } else {
                // the last path separator is first character
                // so this is a file in the root directory
                _dirname.push_back(PATH_SEPARATOR);
            }
        } else {
            // relative file path in current directory
            _basename.assign(_path);
        }

        if(!is_standard_stream()) {
            // split extension
            position = _basename.find_last_of(EXTENSION_SEPARATOR);
            if(position != string::npos) {
                if(position + 1 < _basename.size()) {
                    _extension.assign(_basename, position + 1, string::npos);
                }
                _basename.resize(position);

                // if there is a second extension
                // and the first is a compression marker
                position = _basename.find_last_of(EXTENSION_SEPARATOR);
                if(position != string::npos && (_extension == "gz" || _extension == "bz2" || _extension == "xz")) {
                    _compression.assign(_extension);
                    _extension.clear();
                    if(position + 1 < _basename.size()) {
                        _extension.assign(_basename, position + 1, string::npos);
                    }
                    _basename.resize(position);
                }
            }

            if(!_extension.empty() && _type == FormatType::UNKNOWN) {
                from_string(_extension, _type);
            }
        }
        refresh();
    }
};
void URL::set_basename(const string& name) {
    _basename.assign(name);
    refresh();
};
void URL::set_dirname(const string& directory) {
    _dirname.assign(directory);
    expand(_dirname);
    refresh();
};
void URL::set_compression(const string& compression) {
    _compression.assign(compression);
    refresh();
};
void URL::set_type(const string& type) {
    from_string(type, _type);
    if(!is_standard_stream()) {
        if(_type != FormatType::UNKNOWN && _extension.empty()) {
            decode_extension(_type);
        }
        refresh();
    }
};
void URL::set_type(const FormatType type, const bool force) {
    _type = type;
    if(!is_standard_stream()) {
        if(_type != FormatType::UNKNOWN && (force || _extension.empty())) {
            decode_extension(_type);
        }
        refresh();
    }
};
void URL::relocate(const URL& base) {
    if(!base._dirname.empty() && !is_absolute()) {
        string joined;
        joined.append(base._dirname);
        if(!_dirname.empty()) {
            if(joined.back() != PATH_SEPARATOR) {
                joined.push_back(PATH_SEPARATOR);
            }
            joined.append(_dirname);
        }
        _dirname.assign(joined);
        refresh();
    }
};
bool URL::is_readable() const {
    if(is_stdin()) {
        return true;
    } else if(is_stdout() || is_stderr() || is_null()) {
        return false;
    } else {
        return access(_path.c_str(), R_OK) != -1;
    }
};
bool URL::is_writable() const {
    if(is_stdin()) {
        return false;
    } else if(is_stdout() || is_stderr() || is_null()) {
        return true;
    } else {
        return access(_path.c_str(), F_OK) != -1 ? access(_path.c_str(), W_OK) != -1 : access(_dirname.c_str(), W_OK) != -1;
    }
};
const char* const URL::c_str() const {
    return _path.c_str();
};
const size_t URL::size() const {
    return _path.size();
};
bool URL::operator==(const URL &other) const {
    return _path == other._path;
};
URL& URL::operator=(const URL& other) {
    if(this != &other) {
        clear();
        _path.assign(other._path);
        _basename.assign(other._basename);
        _dirname.assign(other._dirname);
        _extension.assign(other._extension);
        _compression.assign(other._compression);
        _type = other._type;
    }
    return *this;
};
URL::operator string() const {
    return string(_path);
};
void URL::refresh() {
    _path.clear();
    if(!_dirname.empty()) {
        _path.append(_dirname);
    }
    if(!_basename.empty()) {
        if(!_path.empty() && _path.back() != PATH_SEPARATOR) {
            _path.push_back(PATH_SEPARATOR);
        }
        _path.append(_basename);
    }
    if(!is_standard_stream()) {
        if(!_extension.empty()) {
            _path.push_back(EXTENSION_SEPARATOR);
            _path.append(_extension);
        }
        if(!_compression.empty()) {
            _path.push_back(EXTENSION_SEPARATOR);
            _path.append(_compression);
        }
    }
};
void URL::decode_extension(const FormatType& type) {
    _extension.clear();
    to_string(type, _extension);
};
void URL::expand(string& path) {
    if(!path.empty()) {
        string resolved;
        string variable;
        char* value(NULL);
        uint64_t position(0);
        while(position < path.size()) {
            const char& c = path[position];
            switch(c) {
                case '~':
                    if(resolved.empty()) {
                        value = getenv("HOME");
                        if(value != NULL) {
                            resolved.append(value);
                        } else {
                            resolved.push_back(c);
                        }
                    } else {
                        resolved.push_back(c);
                    }
                    break;
                case '$':
                    if(variable.empty()) {
                        variable.push_back(c);
                    } else {
                        resolved.push_back(c);
                    }
                    break;
                case '{':
                    if(variable.size() == 1) {
                        variable.push_back(c);
                    } else {
                        resolved.push_back(c);
                    }
                    break;
                case '}':
                    if(variable.empty()) {
                        resolved.push_back(c);
                    } else if(variable.size() == 1) {
                        resolved.append(variable);
                        resolved.push_back(c);
                        variable.clear();
                    } else {
                        value = getenv(variable.c_str() + 2);
                        if(value != NULL) {
                            resolved.append(value);
                        }
                        variable.clear();
                    }
                    break;
                case '\0':
                    throw InternalError("unexpected string terminal encountered in " + string(*this));
                    break;
                default:
                    if(variable.empty()) {
                        resolved.push_back(c);
                    } else if(variable.size() == 1) {
                        resolved.append(variable);
                        resolved.push_back(c);
                        variable.clear();
                    } else {
                        variable.push_back(c);
                    }
            }
            position++;
        }
        if(variable.empty()) {
        } else if(variable.size() == 1) {
            resolved.append(variable);
            variable.clear();
        } else {
            string message("unterminated environment variable in path ");
            message += path;
            message += " at position ";
            message += to_string(position - variable.size());
            throw ConfigurationError (message);
        }
        path.assign(resolved);
    }
};
bool operator<(const URL& lhs, const URL& rhs) {
    return lhs._path < rhs._path;
};
ostream& operator<<(ostream& o, const URL& url) {
    o << url._path;
    return o;
};
bool decode_directory_url_by_key(const Value::Ch* key, URL& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsString()) {
            value.parse_directory(string(element->value.GetString(), element->value.GetStringLength()));
            return true;
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
bool decode_file_url_by_key(const Value::Ch* key, URL& value, const IoDirection& direction, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        value.clear();
        if(element->value.IsString() || element->value.IsObject()) {
            try {
                if(element->value.IsString()) {
                    value.parse_file(string(element->value.GetString(), element->value.GetStringLength()), direction);
                    return true;
                } else {
                    string buffer;
                    if(decode_value_by_key< string >("path", buffer, element->value)) {
                        value.parse_file(buffer, direction);
                        if(decode_value_by_key< string >("type", buffer, element->value)) {
                            value.set_type(buffer);
                        }
                        if(decode_value_by_key< string >("compression", buffer, element->value)) {
                            value.set_compression(buffer);
                        }
                        return true;
                    } else { throw ConfigurationError("URL element must contain a non empty path element"); }
                }
            } catch(ConfigurationError& e) {
                value.clear();
                throw e;
            }
        } else { throw ConfigurationError("URL element must be either a string or a dictionary"); }
    }
    return false;
};
bool decode_file_url_list_by_key(const Value::Ch* key, list< URL >& value, const Value& container, const IoDirection& direction) {
    Value::ConstMemberIterator collection = container.FindMember(key);
    if(collection != container.MemberEnd()) {
        if(!collection->value.IsNull()) {
            if(collection->value.IsArray() && !collection->value.Empty()) {
                value.clear();
                for(const auto& element : collection->value.GetArray()) {
                    if(element.IsString() || element.IsObject()) {
                        try {
                            if(element.IsString()) {
                                value.emplace_back(string(element.GetString(), element.GetStringLength()), direction);
                            } else {
                                string buffer;
                                if(decode_value_by_key< string >("path", buffer, element)) {
                                    URL url(buffer, direction);
                                    if(decode_value_by_key< string >("type", buffer, element)) {
                                        url.set_type(buffer);
                                    }
                                    if(decode_value_by_key< string >("compression", buffer, element)) {
                                        url.set_compression(buffer);
                                    }
                                    value.emplace_back(url);
                                } else { throw ConfigurationError("URL element must contain a non empty path element"); }
                            }
                        } catch(ConfigurationError& e) {
                            value.clear();
                            throw e;
                        }
                    } else { throw ConfigurationError("URL element must be either a string or a dictionary"); }
                }
                return true;
            }
        }
    }
    return false;
};
bool encode_key_value(const string& key, const URL& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value v(value.c_str(), value.size(), document.GetAllocator());
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
bool encode_key_value(const string& key, const list< URL >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v.c_str(), v.size(), document.GetAllocator()).Move(), document.GetAllocator());
        }
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(k.Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
void encode_element(const URL& value, Value& container, Document& document) {
    Value v(value.c_str(), value.size(), document.GetAllocator());
    container.PushBack(v.Move(), document.GetAllocator());
};
