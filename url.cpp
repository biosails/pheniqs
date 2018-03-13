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

#include "url.h"

ostream& operator<<(ostream& o, const IoDirection& direction) {
    switch (direction) {
        case IoDirection::IN:   o << "in";  break;
        case IoDirection::OUT:  o << "out"; break;
    }
    return o;
};

URL::URL() :
    _type(FormatType::UNKNOWN) {
};
URL::URL(const URL& other) :
    _path(other._path),
    _name(other._name),
    _directory(other._directory),
    _extension(other._extension),
    _compression(other._compression),
    _type(other._type) {
};
URL::URL(const string& path, const IoDirection& direction) : 
    _type(FormatType::UNKNOWN) {
    parse(path, direction);
};
void URL::parse(const string& path, const IoDirection& direction) {
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
            }
        } else {
            if(_path == "/dev/stdin" || 
                _path == "/dev/fd/0" || 
                _path == "/proc/self/fd/0") {
                _path.assign(CANONICAL_STDIN_PATH);
            } else if(_path == "/dev/stdout" ||
                _path == "/dev/fd/1" ||
                _path == "/proc/self/fd/1") {
                _path.assign(CANONICAL_STDOUT_PATH);
            } else if(_path == "/dev/stderr" ||
                _path == "/dev/fd/2" ||
                _path == "/proc/self/fd/2") {
                _path.assign(CANONICAL_STDERR_PATH);
            }
        }

        // split the path into dirname and basename
        auto position = _path.find_last_of(PATH_SEPARATOR);
        if(position != string::npos) {
            // the path separator is present
            if(position  + 1 < _path.size()) {
                _name.assign(_path, position + 1, string::npos);
            }
            if(position > 0) {
                // the last path separator is not the first character
                _directory.assign(_path, 0, position);
            } else {
                // the last path separator is first character
                // so this is a file in the root directory
                _directory.push_back(PATH_SEPARATOR);
            }
        } else {
            // relative file path in current directory
            _name.assign(_path);
        }

        if(!is_standard_stream()) {
            // split extension
            position = _name.find_last_of(EXTENSION_SEPARATOR);
            if(position != string::npos) {
                if(position + 1 < _name.size()) {
                    _extension.assign(_name, position + 1, string::npos);
                }
                _name.resize(position);

                // if there is a second extension 
                // and the first is a compression marker
                position = _name.find_last_of(EXTENSION_SEPARATOR);
                if(position != string::npos && (_extension == "gz" || _extension == "bz2" || _extension == "xz")) {
                    _compression.assign(_extension);
                    _extension.clear();
                    if(position + 1 < _name.size()) {
                        _extension.assign(_name, position + 1, string::npos);
                    }
                    _name.resize(position);
                }
            }

            if(!_extension.empty() && _type == FormatType::UNKNOWN) {
                _extension >> _type;
            }
        }
        refresh();
    }
};
void URL::set_name(const string& name) {
    _name.assign(name);
    refresh();
};
void URL::set_directory(const string& directory) {
    _directory.assign(directory);
    expand(_directory);
    refresh();
};
void URL::set_compression(const string& compression) {
    _compression.assign(compression);
    refresh();
};
void URL::set_type(const string& type) {
    type >> _type;
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
    if(!base._directory.empty() && !is_absolute()) {
        string joined;
        joined.append(base._directory);
        if(!_directory.empty()) {
            if(joined.back() != PATH_SEPARATOR) {
                joined.push_back(PATH_SEPARATOR);
            }
            joined.append(_directory);
        }
        _directory.assign(joined);
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
        return access(_path.c_str(), F_OK) != -1 ? access(_path.c_str(), W_OK) != -1 : access(_directory.c_str(), W_OK) != -1;
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
        _name.assign(other._name);
        _directory.assign(other._directory);
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
    if(!_directory.empty()) {
        _path.append(_directory);
    }
    if(!_name.empty()) {
        if(!_path.empty() && _path.back() != PATH_SEPARATOR) {
            _path.push_back(PATH_SEPARATOR);
        }
        _path.append(_name);
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
    _extension << type;
};
void URL::expand(string& path) {
    if(!path.empty()) {
        string resolved;
        string variable;
        char* value = NULL;
        size_t position = 0;
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
void decode_directory_by_key(const Value::Ch* key, URL& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsString()) {
            value.set_directory(string(element->value.GetString(), element->value.GetStringLength()));
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
};
void decode_url_by_key(const Value::Ch* key, URL& value, const IoDirection& direction, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        value.clear();
        if(element->value.IsString() || element->value.IsObject()) {
            try {
                if(element->value.IsString()) {
                    value.parse(string(element->value.GetString(), element->value.GetStringLength()), direction);

                } else {
                    string buffer;
                    decode_string_by_key("path", buffer, element->value);
                    if(!buffer.empty()) {
                        value.parse(buffer, direction);
                    } else { throw ConfigurationError("URL element must contain a non empty path element"); }

                    buffer.clear();
                    decode_string_by_key("type", buffer, element->value);
                    if(!buffer.empty()) {
                        value.set_type(buffer);
                    }

                    buffer.clear();
                    decode_string_by_key("compression", buffer, element->value);
                    if(!buffer.empty()) {
                        value.set_compression(buffer);
                    }
                }
            } catch(ConfigurationError& e) {
                value.clear();
                throw e;
            }
        } else { throw ConfigurationError("URL element must be either a string or a dictionary"); }
    }
};
void encode_key_value(const string& key, const URL& value, Value& node, Document& document) {
    if(!value.empty()) {
        Value v(value.c_str(), value.size(), document.GetAllocator());
        Value k(key.c_str(), key.size(), document.GetAllocator());
        node.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};