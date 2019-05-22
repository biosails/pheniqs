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

string to_string(const IoDirection& value) {
    string result;
    switch (value) {
        case IoDirection::IN:   result.assign("in");    break;
        case IoDirection::OUT:  result.assign("out");   break;
        default:                                        break;
    }
    return result;
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
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const IoDirection& value, Value& container, Document& document) {
    encode_key_value(key, to_string(value), container, document);
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

string to_string(const FormatType& value) {
    string result;
    switch (value) {
        case FormatType::NONE:      result.assign("none");   break;
        case FormatType::FASTQ:     result.assign("fastq");  break;
        case FormatType::SAM:       result.assign("sam");    break;
        case FormatType::BAM:       result.assign("bam");    break;
        case FormatType::BAI:       result.assign("bai");    break;
        case FormatType::CRAM:      result.assign("cram");   break;
        case FormatType::CRAI:      result.assign("crai");   break;
        case FormatType::VCF:       result.assign("vcf");    break;
        case FormatType::BCF:       result.assign("bcf");    break;
        case FormatType::CSI:       result.assign("csi");    break;
        case FormatType::GZI:       result.assign("gzi");    break;
        case FormatType::TBI:       result.assign("tbi");    break;
        case FormatType::BED:       result.assign("bed");    break;
        case FormatType::JSON:      result.assign("json");   break;
        default:                                             break;
    }
    return result;
};
bool from_string(const char* value, FormatType& result) {
         if(value == NULL)              result = FormatType::UNKNOWN;
    else if(!strcmp(value, "none"))     result = FormatType::NONE;
    else if(!strcmp(value, "fastq"))    result = FormatType::FASTQ;
    else if(!strcmp(value, "fq"))       result = FormatType::FASTQ;
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
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const FormatType& value, Value& container, Document& document) {
    encode_key_value(key, to_string(value), container, document);
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

string to_string(const FormatCompression& value) {
    string result;
    switch (value) {
        case FormatCompression::NONE:      result.assign("none");   break;
        case FormatCompression::GZIP:      result.assign("gz");     break;
        case FormatCompression::BGZF:      result.assign("bgzf");   break;
        case FormatCompression::BZ2:       result.assign("bz2");    break;
        case FormatCompression::XZ:        result.assign("xz");     break;
        default:                                                    break;
    }
    return result;
};
bool from_string(const char* value, FormatCompression& result) {
         if(value == NULL)              result = FormatCompression::UNKNOWN;
    else if(!strcmp(value, "none"))     result = FormatCompression::NONE;
    else if(!strcmp(value, "gz"))       result = FormatCompression::GZIP;
    else if(!strcmp(value, "bgzf"))     result = FormatCompression::BGZF;
    else if(!strcmp(value, "bz2"))      result = FormatCompression::BZ2;
    else if(!strcmp(value, "xz"))       result = FormatCompression::XZ;
    else                                result = FormatCompression::UNKNOWN;
    return (result == FormatCompression::UNKNOWN ? false : true);
};
bool from_string(const string& value, FormatCompression& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const FormatCompression& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const FormatCompression& value, Value& container, Document& document) {
    encode_key_value(key, to_string(value), container, document);
};
template<> bool decode_value_by_key< FormatCompression >(const Value::Ch* key, FormatCompression& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> FormatCompression decode_value_by_key(const Value::Ch* key, const Value& container) {
    FormatCompression value(FormatCompression::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

string to_string(const CompressionLevel& value) {
    string result;
    switch (value) {
    case CompressionLevel::LEVEL_0: result.assign("0"); break;
    case CompressionLevel::LEVEL_1: result.assign("1"); break;
    case CompressionLevel::LEVEL_2: result.assign("2"); break;
    case CompressionLevel::LEVEL_3: result.assign("3"); break;
    case CompressionLevel::LEVEL_4: result.assign("4"); break;
    case CompressionLevel::LEVEL_5: result.assign("5"); break;
    case CompressionLevel::LEVEL_6: result.assign("6"); break;
    case CompressionLevel::LEVEL_7: result.assign("7"); break;
    case CompressionLevel::LEVEL_8: result.assign("8"); break;
    case CompressionLevel::LEVEL_9: result.assign("9"); break;
    default:                                                break;
    }
    return result;
};
bool from_string(const char* value, CompressionLevel& result) {
         if(value == NULL)          result = CompressionLevel::UNKNOWN;
    else if(!strcmp(value, "0"))    result = CompressionLevel::LEVEL_0;
    else if(!strcmp(value, "1"))    result = CompressionLevel::LEVEL_1;
    else if(!strcmp(value, "2"))    result = CompressionLevel::LEVEL_2;
    else if(!strcmp(value, "3"))    result = CompressionLevel::LEVEL_3;
    else if(!strcmp(value, "4"))    result = CompressionLevel::LEVEL_4;
    else if(!strcmp(value, "5"))    result = CompressionLevel::LEVEL_5;
    else if(!strcmp(value, "6"))    result = CompressionLevel::LEVEL_6;
    else if(!strcmp(value, "7"))    result = CompressionLevel::LEVEL_7;
    else if(!strcmp(value, "8"))    result = CompressionLevel::LEVEL_8;
    else if(!strcmp(value, "9"))    result = CompressionLevel::LEVEL_9;
    else                            result = CompressionLevel::UNKNOWN;
    return (result == CompressionLevel::UNKNOWN ? false : true);
};
bool from_string(const string& value, CompressionLevel& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const CompressionLevel& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const CompressionLevel& value, Value& container, Document& document) {
    encode_key_value(key, to_string(value), container, document);
};
template<> bool decode_value_by_key< CompressionLevel >(const Value::Ch* key, CompressionLevel& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> CompressionLevel decode_value_by_key(const Value::Ch* key, const Value& container) {
    CompressionLevel value(CompressionLevel::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

string to_string(const URLQueryParameter& value) {
    string result;
    switch (value) {
    case URLQueryParameter::FORMAT_TYPE:        result.assign("format");        break;
    case URLQueryParameter::FORMAT_COMPRESSION: result.assign("compression");   break;
    case URLQueryParameter::COMPRESSION_LEVEL:  result.assign("level");         break;
    default:                                                                    break;
    }
    return result;
};
bool from_string(const char* value, URLQueryParameter& result) {
         if(value == NULL)                  result = URLQueryParameter::UNKNOWN;
    else if(!strcmp(value, "format"))       result = URLQueryParameter::FORMAT_TYPE;
    else if(!strcmp(value, "compression"))  result = URLQueryParameter::FORMAT_COMPRESSION;
    else if(!strcmp(value, "level"))        result = URLQueryParameter::COMPRESSION_LEVEL;
    else                                    result = URLQueryParameter::UNKNOWN;
    return (result == URLQueryParameter::UNKNOWN ? false : true);
};
bool from_string(const string& value, URLQueryParameter& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const URLQueryParameter& value) {
    o << to_string(value);
    return o;
};

URL::URL() :
    _format_type(FormatType::UNKNOWN),
    _format_compression(FormatCompression::UNKNOWN),
    _compression_level(CompressionLevel::UNKNOWN) {
};
URL::URL(const URL& other) :
    _encoded(other._encoded),
    _path(other._path),
    _basename(other._basename),
    _dirname(other._dirname),
    _query(other._query),
    _format_type(other._format_type),
    _format_compression(other._format_compression),
    _compression_level(other._compression_level) {
};
URL::URL(const string& encoded) :
    _format_type(FormatType::UNKNOWN),
    _format_compression(FormatCompression::UNKNOWN),
    _compression_level(CompressionLevel::UNKNOWN) {
    parse(encoded);
};
void URL::parse(const string& encoded) {
    clear();
    if(!encoded.empty()) {
        _encoded.assign(encoded);

        /* split the path into dirname and basename using the last instance of the path separator */
        auto position = _encoded.find_last_of('/');
        if(position != string::npos) {
            /* the path separator is present */
            if(position  + 1 < _encoded.size()) {
                /* if the path separator is not the last character the basename is not empty */
                _basename.assign(_encoded, position + 1, string::npos);
            }
            if(position > 0) {
                /* if the last path separator is not the first character dirname is not the root */
                _dirname.assign(_encoded, 0, position);
            } else {
                /* if the last path separator is the first character dirname is the root */
                _dirname.push_back('/');
            }
        } else {
            /* if the path separator is not present this is a relative basename in the current directory */
            _basename.assign(_encoded);
        }

        /* split the query from the basename if the query separator is present */
        position = _basename.find_first_of('?');
        if(position != string::npos) {
            /* if the query separator is not the last character the query is not empty */
            if(position  + 1 < _basename.size()) {
                _query.assign(_basename, position + 1, string::npos);
            }
            _basename.erase(position);
        }
        if(_basename == "." || _basename == "..") {
            /* basename is a reference to the current or parent directory */
            if(!_dirname.empty()) { _dirname.push_back('/'); }
            _dirname.append(_basename);
            _basename.erase(0);
        } else {
            /* if the extension separator is present try and extract an implicit compression and type */
            string buffer(_basename);
            string extension;
            position = buffer.find_last_of('.');
            if(position != string::npos) {
                if(position > 0) {
                    if(position + 2 < buffer.size()) {
                        extension.assign(buffer, position + 1, string::npos);
                        buffer.erase(position);
                        /* if the trailing extension is a compression extension
                           and there is another trailing extension implicitly set the compression
                           and decode the trailing extension */
                        if(extension == "gz" || extension == "bz2" || extension == "xz") {
                            /* set the implicit compression */
                            from_string(extension, _format_compression);
                            position = buffer.find_last_of('.');
                            if(position != string::npos) {
                                extension.clear();
                                if(position + 2 < buffer.size()) {
                                    extension.assign(buffer, position + 1, string::npos);
                                }
                            }
                        }
                        if(!extension.empty()) {
                            /* set the implicit format type */
                            from_string(extension, _format_type);
                        }
                    }
                }
            }
        }
        parse_query();
        refresh();
    }
};
void URL::parse_query() {
    if(!_query.empty()) {
        string key;
        string value;
        uint8_t segment = 0;
        const char* buffer = _query.c_str();
        for(string::size_type position = 0; position < _query.size(); ++position) {
            char c(*(buffer + position));
            switch(c) {
                case '=':
                    segment = 1;
                    break;
                case '&':
                    apply_query_parameter(key, value);
                    key.clear();
                    value.clear();
                    segment = 0;
                    break;
                default:
                    switch(segment) {
                        case 0:
                            key.push_back(c);
                            break;
                        case 1:
                            value.push_back(c);
                            break;
                        default:
                            break;
                    };
                    break;
            };
        }
        apply_query_parameter(key, value);
    }
};
void URL::apply_query_parameter(const string& key, const string& value) {
    if(!key.empty() && !value.empty()) {
        URLQueryParameter name;
        from_string(key, name);
        switch(name) {
            case URLQueryParameter::FORMAT_COMPRESSION: {
                from_string(value, _format_compression);
                break;
            };
            case URLQueryParameter::FORMAT_TYPE: {
                from_string(value, _format_type);
                break;
            };
            case URLQueryParameter::COMPRESSION_LEVEL: {
                from_string(value, _compression_level);
                break;
            };
            case URLQueryParameter::UNKNOWN: {
                /* silently ignore unknown parameters */
                break;
            };
        }
    }
};
void URL::refresh() {
    /* rebuild path */
    _path.clear();
    if(!_dirname.empty()) {
        _path.append(_dirname);
    }
    if(!_basename.empty()) {
        if(!_path.empty() && _path.back() != '/') {
            _path.push_back('/');
        }
        _path.append(_basename);
    }

    /* rebuild query */
    _query.clear();
    if(_format_type != FormatType::UNKNOWN) {
        if (!_query.empty()) { _query.push_back('&'); }
        _query.append(to_string(URLQueryParameter::FORMAT_TYPE));
        _query.push_back('=');
        _query.append(to_string(_format_type));

        /* make sure only valid compression is specified for each format */
        switch(_format_type) {
            case FormatType::SAM: {
                _format_compression = FormatCompression::NONE;
                _compression_level = CompressionLevel::UNKNOWN;
                break;
            };
            case FormatType::BAM:
            case FormatType::FASTQ: {
                switch(_format_compression) {
                    case FormatCompression::NONE: {
                        _compression_level = CompressionLevel::UNKNOWN;
                        break;
                    };
                    case FormatCompression::GZIP: {
                        break;
                    };
                    case FormatCompression::BGZF: {
                        break;
                    };
                    default:
                        _format_compression = FormatCompression::UNKNOWN;
                        break;
                }
                break;
            };
            case FormatType::CRAM: {
                _format_compression = FormatCompression::UNKNOWN;
                break;
            };
            case FormatType::JSON: {
                _format_compression = FormatCompression::NONE;
                _compression_level = CompressionLevel::UNKNOWN;
                break;
            };
            default: {
                _format_compression = FormatCompression::UNKNOWN;
                break;
            };
        }

        if(_format_compression != FormatCompression::UNKNOWN) {
            if (!_query.empty()) { _query.push_back('&'); }
            _query.append(to_string(URLQueryParameter::FORMAT_COMPRESSION));
            _query.push_back('=');
            _query.append(to_string(_format_compression));
        }

        if(_compression_level != CompressionLevel::UNKNOWN) {
            if (!_query.empty()) { _query.push_back('&'); }
            _query.append(to_string(URLQueryParameter::COMPRESSION_LEVEL));
            _query.push_back('=');
            _query.append(to_string(_compression_level));
        }
    }

    /* rebuild encoded */
    _encoded.clear();
    if(!_path.empty()) {
        _encoded.append(_path);
    }
    if(!_query.empty()) {
        _encoded.push_back('?');
        _encoded.append(_query);
    }
};
void URL::set_type(const FormatType type) {
    _format_type = type;
    refresh();
};
void URL::set_compression(const FormatCompression& compression) {
    _format_compression = compression;
    refresh();
};
void URL::set_compression_level(const CompressionLevel& level) {
    _compression_level = level;
    refresh();
};
void URL::override_query(const URL& other) {
    if(other._format_type != FormatType::UNKNOWN) {
        _format_type = other._format_type;
    }
    if(other._format_compression != FormatCompression::UNKNOWN) {
        _format_compression = other._format_compression;
    }
    if(other._compression_level != CompressionLevel::UNKNOWN) {
        _compression_level = other._compression_level;
    }
    refresh();
};
void URL::relocate_child(const URL& base) {
    if(!base._path.empty() && !is_absolute()) {
        string joined;
        joined.append(base._path);
        if(!_dirname.empty()) {
            if(joined.back() != '/') {
                joined.push_back('/');
            }
            joined.append(_dirname);
        }
        _dirname.assign(joined);
        refresh();
    }
};
void URL::relocate_sibling(const URL& base) {
    if(!base._dirname.empty() && !is_absolute()) {
        string joined;
        joined.append(base._dirname);
        if(!_dirname.empty()) {
            if(joined.back() != '/') {
                joined.push_back('/');
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
    } else if(is_stdout() || is_stderr() || is_dev_null()) {
        return false;
    } else {
        return !access(_path.c_str(), R_OK);
    }
};
bool URL::is_writable() const {
    if(is_stdin()) {
        return false;
    } else if(is_stdout() || is_stderr() || is_dev_null()) {
        return true;
    } else {
        return !access(_path.c_str(), F_OK) ? !access(_path.c_str(), W_OK) : !access(_dirname.c_str(), W_OK);
    }
};
bool URL::operator==(const URL &other) const {
    return _encoded == other._encoded;
};
URL& URL::operator=(const URL& other) {
    if(this != &other) {
        clear();
        _encoded.assign(other._encoded);
        _path.assign(other._path);
        _basename.assign(other._basename);
        _dirname.assign(other._dirname);
        _query.assign(other._query);
        _format_type = other._format_type;
        _format_compression = other._format_compression;
        _compression_level = other._compression_level;
    }
    return *this;
};
URL::operator string() const {
    return string(_encoded);
};
string URL::description() const {
    string description;
    description.append("URL                : ");
    description.append(_encoded);
    description.append("\n");

    description.append("Path               : ");
    description.append(_path);
    description.append("\n");

    description.append("Basename           : ");
    description.append(_basename);
    description.append("\n");

    description.append("Dirname            : ");
    description.append(_dirname);
    description.append("\n");

    description.append("Query              : ");
    description.append(_query);
    description.append("\n");

    description.append("Type               : ");
    description.append(to_string(_format_type));
    description.append("\n");

    description.append("Compression        : ");
    description.append(to_string(_format_compression));
    description.append("\n");

    description.append("Compression level  : ");
    description.append(to_string(_compression_level));
    description.append("\n");

    return description;
};

bool operator<(const URL& lhs, const URL& rhs) {
    return lhs._encoded < rhs._encoded;
};
ostream& operator<<(ostream& o, const URL& url) {
    o << url._encoded;
    return o;
};

string& expand_shell(string& expression) {
    if(!expression.empty()) {
        string resolved;
        string variable;
        char* value(NULL);
        size_t position(0);
        while(position < expression.size()) {
            const char& c = expression[position];
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
                    throw InternalError("unexpected string terminator encountered in " + expression);
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
            ++position;
        }
        if(!variable.empty()) {
            if(variable.size() == 1) {
                resolved.append(variable);
                variable.clear();
            } else {
                string message("unterminated environment variable in expression ");
                message += expression;
                message += " at position ";
                message += to_string(position - variable.size());
                throw ConfigurationError (message);
            }
        }
        expression.assign(resolved);
    }
    return expression;
};
void normalize_standard_stream(string& path, const IoDirection& direction) {
    if(path == STANDARD_STREAM_ALIAS) {
        switch(direction) {
            case IoDirection::IN: {
                path.assign(CANONICAL_STDIN_PATH);
                break;
            };
            case IoDirection::OUT: {
                path.assign(CANONICAL_STDOUT_PATH);
                break;
            };
            case IoDirection::UNKNOWN: {
                throw ConfigurationError("can not interpret standard stream alias " + path + " without a declared IO direction");
                break;
            };
        }
    } else {
        if(path == "/dev/stdin" || path == "/dev/fd/0" || path == "/proc/self/fd/0") {
            if(direction == IoDirection::OUT) {
                throw ConfigurationError("can not use " + path + " for output");
            }
            path.assign(CANONICAL_STDIN_PATH);

        } else if(path == "/dev/stdout" || path == "/dev/fd/1" || path == "/proc/self/fd/1") {
            if(direction == IoDirection::IN) {
                throw ConfigurationError("can not use " + path + " for input");
            }
            path.assign(CANONICAL_STDOUT_PATH);

        } else if(path == "/dev/stderr" || path == "/dev/fd/2" || path == "/proc/self/fd/2") {
            if(direction == IoDirection::IN) {
                throw ConfigurationError("can not use " + path + " for input");
            }
            path.assign(CANONICAL_STDERR_PATH);
        }
    }
};
void standardize_url_value(Value& container, Document& document, const IoDirection& direction) {
    if(!container.IsNull()) {
        if(container.IsString()) {
            string buffer(container.GetString(), container.GetStringLength());
            expand_shell(buffer);
            URL url(buffer);
            buffer.assign(url.path());
            normalize_standard_stream(buffer, direction);
            if(!url.query().empty()) {
                buffer.push_back('?');
                buffer.append(url.query());
            }
            url.parse(buffer);
            container.SetString(url.encoded().c_str(), url.encoded().size(), document.GetAllocator());
        }
    }
};
void standardize_url_value_by_key(const Value::Ch* key, Value& container, Document& document, const IoDirection& direction) {
    if(container.IsObject()) {
        Value::MemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            standardize_url_value(reference->value, document, direction);
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
void standardize_url_array_by_key(const Value::Ch* key, Value& container, Document& document, const IoDirection& direction) {
    if(container.IsObject()) {
        Value::MemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(reference->value.IsArray()) {
                for(auto& record : reference->value.GetArray()) {
                    standardize_url_value(record, document, direction);
                }
            } else { throw ConfigurationError(string(key) + " is not an array");  }
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
void relocate_url_array_by_key(const Value::Ch* key, Value& container, Document& document, const URL& base) {
    list< URL > value;
    if(decode_value_by_key< list< URL > >(key, value, container)) {
        for(auto& url : value) {
            url.relocate_child(base);
        }
    }
    encode_key_value(key, value, container, document);
};

template<> URL decode_value(const Value& container) {
    if(!container.IsNull()) {
        if(container.IsString()) {
            URL value;
            string buffer;
            buffer.assign(container.GetString(), container.GetStringLength());
            value.parse(buffer);
            return value;
        } else { throw ConfigurationError("URL element must be a string"); }
    } else { throw ConfigurationError("URL element is null"); }
};
template<> URL decode_value_by_key(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            return decode_value< URL >(reference->value);
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
template<> bool decode_value< URL >(URL& value, const Value& container) {
    if(!container.IsNull()) {
        if(container.IsString()) {
            try {
                string buffer(container.GetString(), container.GetStringLength());
                value.parse(buffer);
                return true;
            } catch(ConfigurationError& e) {
                value.clear();
                throw e;
            }
        } else { throw ConfigurationError("URL element must be a string"); }
    }
    return false;
};
template<> list< URL > decode_value_by_key(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(!reference->value.IsNull()) {
                list< URL > value;
                if(reference->value.IsArray()) {
                    for(const auto& element : reference->value.GetArray()) {
                        value.emplace_back(decode_value< URL >(element));
                    }
                } else if(reference->value.IsString()) {
                    value.emplace_back(decode_value< URL >(reference->value));
                }
                return value;
            } else { throw ConfigurationError(string(key) + " is null"); }
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
template<> bool decode_value_by_key< URL >(const Value::Ch* key, URL& value, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            return decode_value< URL >(value, reference->value);
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
template<> bool decode_value_by_key< list< URL > >(const Value::Ch* key, list< URL >& value, const Value& container) {
    bool result(false);
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd() && !reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                value.clear();
                for(const auto& element : reference->value.GetArray()) {
                    value.emplace_back(decode_value< URL >(element));
                    result = true;
                }
            } else { throw ConfigurationError(string(key) + " element is not an array"); }
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return result;
};
void encode_value(const URL& value, Value& container, Document& document) {
    container.SetString(value.encoded().c_str(), value.encoded().size(), document.GetAllocator());
};
bool encode_key_value(const string& key, const URL& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value element;
        encode_value(value, element, document);
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
bool encode_key_value(const string& key, const list< URL >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& url : value) {
            Value element;
            encode_value(url, element, document);
            array.PushBack(element.Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
