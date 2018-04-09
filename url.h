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

#ifndef PHENIQS_URL_H
#define PHENIQS_URL_H

#include <stdio.h>
#include <iomanip>
#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "error.h"
#include "json.h"
#include "constant.h"

using std::hash;
using std::setw;
using std::endl;
using std::list;
using std::cerr;
using std::cout;
using std::fixed;
using std::uint64_t;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;
using std::numeric_limits;

#define STANDARD_STREAM_ALIAS "-"
#define CANONICAL_STDIN_PATH "/dev/stdin"
#define CANONICAL_STDOUT_PATH "/dev/stdout"
#define CANONICAL_STDERR_PATH "/dev/stderr"
#define CANONICAL_NULL_DEVICE_PATH "/dev/null"
const char PATH_SEPARATOR = '/';

enum class FormatType : uint8_t {
    UNKNOWN,
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
void to_string(const FormatType& value, string& result);
bool from_string(const char* value, FormatType& result);
bool from_string(const string& value, FormatType& result);
ostream& operator<<(ostream& o, const FormatType& value);
void encode_key_value(const string& key, const FormatType& value, Value& container, Document& document);

enum class IoDirection : uint8_t {
    IN,
    OUT,
    UNKNOWN,
};
void to_string(const IoDirection& value, string& result);
bool from_string(const char* value, IoDirection& result);
bool from_string(const string& value, IoDirection& result);
ostream& operator<<(ostream& o, const IoDirection& direction);
void encode_key_value(const string& key, const IoDirection& value, Value& container, Document& document);

class URL {
friend ostream& operator<<(ostream& o, const URL& url);
friend bool operator<(const URL& lhs, const URL& rhs);

public:
    URL();
    URL(const URL& other);
    URL(const string& path);
    URL(const string& path, const bool& is_directory);
    URL(const string& path, const IoDirection& direction);
    void parse_directory(const string& path);
    void parse_file(const string& path, const IoDirection& direction);
    void set_basename(const string& name);
    void set_dirname(const string& directory);
    void set_compression(const string& compression);
    void set_type(const string& type);
    void set_type(const FormatType type, const bool force = false);
    void relocate(const URL& base);
    inline const string& path() const {
        return _path;
    };
    inline const string& basename() const {
        return _basename;
    };
    inline const string& dirname() const {
        return _dirname;
    };
    inline const string& extension() const {
        return _extension;
    };
    inline const string& compression() const {
        return _compression;
    };
    inline const FormatType& type() const {
        return _type;
    };
    inline const FormatKind kind() const {
        switch(_type) {
            case FormatType::SAM:
            case FormatType::BAM:
            case FormatType::BAI:
            case FormatType::CRAM:
            case FormatType::CRAI:
            case FormatType::VCF:
            case FormatType::BCF:
            case FormatType::CSI:
            case FormatType::GZI:
            case FormatType::TBI:
            case FormatType::BED:
                return FormatKind::HTS;
                break;
            case FormatType::FASTQ:
                return FormatKind::FASTQ;
                break;
            default:
                return FormatKind::UNKNOWN;
                break;
        }
    };
    inline bool empty() const {
        return _path.empty();
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
    inline bool is_null() const {
        return _path == CANONICAL_NULL_DEVICE_PATH;
    };
    inline bool is_standard_stream() const {
        return is_stdin() || is_stdout() || is_stderr() || is_null();
    };
    inline bool is_absolute() const {
        return !_dirname.empty() && _dirname[0] == PATH_SEPARATOR;
    };
    inline void clear() {
        _path.clear();
        _basename.clear();
        _dirname.clear();
        _extension.clear();
        _compression.clear();
    };
    bool is_readable() const;
    bool is_writable() const;
    const char* const c_str() const;
    const uint64_t size() const;
    void describe(ostream& o) const;
    bool operator==(const URL& other) const;
    URL& operator=(const URL& other);
    operator string() const;

private:
    string _path;
    string _basename;
    string _dirname;
    string _extension;
    string _compression;
    FormatType _type;

    void refresh();
    void expand(string& path);
    void decode_extension(const FormatType& type);
};
bool operator<(const URL& lhs, const URL& rhs);

namespace std {
    template <> struct hash< URL > {
        uint64_t operator()(const URL& url) const {
            return hash< string >()(url.path());
        };
    };
};
bool decode_directory_url_by_key(const Value::Ch* key, URL& value, const Value& container);
bool decode_file_url_by_key(const Value::Ch* key, URL& value, const IoDirection& direction, const Value& container);
// bool decode_file_url_vector_by_key(const Value::Ch* key, vector< URL >& value, const Value& container, const IoDirection& direction);
bool decode_file_url_list_by_key(const Value::Ch* key, list< URL >& value, const Value& container, const IoDirection& direction);
bool encode_key_value(const string& key, const URL& value, Value& container, Document& document);
// bool encode_key_value(const string& key, const vector< URL >& value, Value& node, Document& document);
bool encode_key_value(const string& key, const list< URL >& value, Value& container, Document& document);
void encode_element(const URL& value, Value& container, Document& document);

#endif /* PHENIQS_URL_H */