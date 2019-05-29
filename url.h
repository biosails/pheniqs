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

#ifndef PHENIQS_URL_H
#define PHENIQS_URL_H

#include "include.h"
#include "json.h"

#define STANDARD_STREAM_ALIAS "-"
#define CANONICAL_STDIN_PATH "/dev/stdin"
#define CANONICAL_STDOUT_PATH "/dev/stdout"
#define CANONICAL_STDERR_PATH "/dev/stderr"
#define CANONICAL_NULL_DEVICE_PATH "/dev/null"

enum class IoDirection : uint8_t {
    IN,
    OUT,
    UNKNOWN,
};
string to_string(const IoDirection& value);
bool from_string(const char* value, IoDirection& result);
bool from_string(const string& value, IoDirection& result);
ostream& operator<<(ostream& o, const IoDirection& direction);
void encode_key_value(const string& key, const IoDirection& value, Value& container, Document& document);

enum class FormatType : uint8_t {
    UNKNOWN,
    NONE,
    FASTQ,
    SAM,
    BAM,
    BAI,
    CRAM,
    CRAI,
    VCF,
    BCF,
    CSI,
    GZI,
    TBI,
    BED,
    JSON,
};
string to_string(const FormatType& value);
bool from_string(const char* value, FormatType& result);
bool from_string(const string& value, FormatType& result);
ostream& operator<<(ostream& o, const FormatType& value);
void encode_key_value(const string& key, const FormatType& value, Value& container, Document& document);
template<> bool decode_value_by_key< FormatType >(const Value::Ch* key, FormatType& value, const Value& container);
template<> FormatType decode_value_by_key< FormatType >(const Value::Ch* key, const Value& container);

enum class FormatCompression : uint8_t {
    UNKNOWN,
    NONE,
    GZIP,
    BGZF,
    BZ2,
    XZ,
};
string to_string(const FormatCompression& value);
bool from_string(const char* value, FormatCompression& result);
bool from_string(const string& value, FormatCompression& result);
ostream& operator<<(ostream& o, const FormatCompression& value);
void encode_key_value(const string& key, const FormatCompression& value, Value& container, Document& document);
template<> bool decode_value_by_key< FormatCompression >(const Value::Ch* key, FormatCompression& value, const Value& container);
template<> FormatCompression decode_value_by_key< FormatCompression >(const Value::Ch* key, const Value& container);

enum class CompressionLevel : uint8_t {
    LEVEL_0,
    LEVEL_1,
    LEVEL_2,
    LEVEL_3,
    LEVEL_4,
    LEVEL_5,
    LEVEL_6,
    LEVEL_7,
    LEVEL_8,
    LEVEL_9,
    UNKNOWN,
};
string to_string(const CompressionLevel& value);
bool from_string(const char* value, CompressionLevel& result);
bool from_string(const string& value, CompressionLevel& result);
ostream& operator<<(ostream& o, const CompressionLevel& value);
void encode_key_value(const string& key, const CompressionLevel& value, Value& container, Document& document);
template<> bool decode_value_by_key< CompressionLevel >(const Value::Ch* key, CompressionLevel& value, const Value& container);
template<> CompressionLevel decode_value_by_key< CompressionLevel >(const Value::Ch* key, const Value& container);

enum class URLQueryParameter : uint8_t {
    UNKNOWN,
    FORMAT_TYPE,
    FORMAT_COMPRESSION,
    COMPRESSION_LEVEL,
};
string to_string(const URLQueryParameter& value);
bool from_string(const char* value, URLQueryParameter& result);
bool from_string(const string& value, URLQueryParameter& result);
ostream& operator<<(ostream& o, const URLQueryParameter& value);

string& expand_shell(string& expression);
void normalize_standard_stream(string& path, const IoDirection& direction);

class URL {
    friend ostream& operator<<(ostream& o, const URL& url);
    friend bool operator<(const URL& lhs, const URL& rhs);

    public:
        URL();
        URL(const URL& other);
        URL(const string& encoded);
        void parse(const string& encoded);
        void set_type(const FormatType type);
        void set_compression(const FormatCompression& compression);
        void set_compression_level(const CompressionLevel& level);
        void override_query(const URL& other);
        void relocate_child(const URL& base);
        void relocate_sibling(const URL& base);
        inline const string& encoded() const {
            return _encoded;
        };
        inline const string& path() const {
            return _path;
        };
        inline const string& basename() const {
            return _basename;
        };
        inline const string& dirname() const {
            return _dirname;
        };
        inline const string& query() const {
            return _query;
        };
        inline const FormatType& type() const {
            return _format_type;
        };
        inline const FormatCompression& implicit_compression() const {
            return  _implicit_compression;
        };
        inline const FormatCompression& explicit_compression() const {
            return  _explicit_compression;
        };
        inline const FormatCompression& compression() const {
            return  _explicit_compression != FormatCompression::UNKNOWN ? _explicit_compression: _implicit_compression;
        };
        inline const CompressionLevel& compression_level() const {
            return _compression_level;
        };
        inline bool empty() const {
            return _encoded.empty();
        };
        inline bool is_file() const {
            return !_basename.empty();
        };
        inline bool is_directory() const {
            return _basename.empty() && !_dirname.empty();
        };
        inline bool is_stdin() const {
            return _path == CANONICAL_STDIN_PATH;
        };
        inline bool is_stdout() const {
            return _path == CANONICAL_STDOUT_PATH;
        };
        inline bool is_stderr() const {
            return _path == CANONICAL_STDERR_PATH;
        };
        inline bool is_dev_null() const {
            return _path == CANONICAL_NULL_DEVICE_PATH;
        };
        inline bool is_standard_stream() const {
            return is_stdin() || is_stdout() || is_stderr() || is_dev_null();
        };
        inline bool is_absolute() const {
            return !_dirname.empty() && _dirname[0] == '/';
        };
        inline const char* const hfile_name() const {
            if(is_file()) {
                if(is_stdout() || is_stdin()) {
                    return STANDARD_STREAM_ALIAS;

                } else if(is_stderr()) {
                    return NULL;

                } else if(is_dev_null()) {
                    return NULL;

                } else {
                    return _path.c_str();
                }
            } else {
                return NULL;
            }
        };
        inline void clear() {
            _encoded.clear();
            _path.clear();
            _basename.clear();
            _dirname.clear();
            _query.clear();
            _format_type = FormatType::UNKNOWN;
            _explicit_compression = FormatCompression::UNKNOWN;
            _implicit_compression = FormatCompression::UNKNOWN;
            _compression_level = CompressionLevel::UNKNOWN;
        };
        bool is_readable() const;
        bool is_writable() const;
        bool operator==(const URL& other) const;
        URL& operator=(const URL& other);
        operator string() const;
        string description() const;

    private:
        string _encoded;
        string _path;
        string _basename;
        string _dirname;
        string _query;
        FormatType _format_type;
        FormatCompression _implicit_compression;
        FormatCompression _explicit_compression;
        CompressionLevel _compression_level;
        void refresh();
        void parse_query();
        void apply_query_parameter(const string& key, const string& value);
        void append_query_parameter(const URLQueryParameter& name, const string& value);

};
bool operator<(const URL& lhs, const URL& rhs);

namespace std {
    template <> struct hash< URL > {
        size_t operator()(const URL& url) const {
            return hash< string >()(url.path());
        };
    };
};

template<> URL decode_value(const Value& container);
template<> URL decode_value_by_key(const Value::Ch* key, const Value& container);
template<> bool decode_value< URL >(URL& value, const Value& container);

template<> list< URL > decode_value_by_key(const Value::Ch* key, const Value& container);
template<> bool decode_value_by_key< URL >(const Value::Ch* key, URL& value, const Value& container);
template<> bool decode_value_by_key< list< URL > >(const Value::Ch* key, list< URL >& value, const Value& container);

void encode_value(const URL& value, Value& container, Document& document);
bool encode_key_value(const string& key, const URL& value, Value& container, Document& document);
bool encode_key_value(const string& key, const list< URL >& value, Value& container, Document& document);

void standardize_url_value_by_key(const Value::Ch* key, Value& container, Document& document, const IoDirection& direction);
void standardize_url_array_by_key(const Value::Ch* key, Value& container, Document& document, const IoDirection& direction);
void relocate_url_by_key(const Value::Ch* key, Value& container, Document& document, const URL& base);
void relocate_url_array_by_key(const Value::Ch* key, Value& container, Document& document, const URL& base);

#endif /* PHENIQS_URL_H */
