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

#include "model.h"

#include <zlib.h>

static inline char* copy_until_tag_end(char* source, const char* end, kstring_t* target) {
    char* position = source;
    while ((*position != '\t' && *position != LINE_BREAK) && position < end) {
        position++;
    }
    uint32_t length = position - source;
    ks_clear(*target);
    if (length > 0) {
        kputsn(source, length, target);
    }
    return position;
};
static inline char* copy_until_linebreak(char* source, const char* end, kstring_t* target) {
    char* position = source;
    while (*position != LINE_BREAK && position < end) {
        position++;
    }
    uint32_t length = position - source;
    ks_clear(*target);
    if (length > 0) {
        kputsn(source, length, target);
    }
    return position;
};
static inline char* skip_to_tab(char* source, const char* end) {
    char* position = source;
    while (*position != '\t' && position < end) {
        position++;
    }
    return position;
};

/* URL
*/
/*
URL::URL() :
    _path({ 0, 0, NULL }),
    _name({ 0, 0, NULL }),
    _directory({ 0, 0, NULL }),
    _extension({ 0, 0, NULL }),
    _compression({ 0, 0, NULL }),
    _type(FormatType::UNKNOWN) {

    initialize();
};
URL::~URL() {
    ks_free(_path);
    ks_free(_name);
    ks_free(_directory);
    ks_free(_extension);
    ks_free(_compression);
};
URL::URL(const URL& other) :
    _path({ 0, 0, NULL }),
    _name({ 0, 0, NULL }),
    _directory({ 0, 0, NULL }),
    _extension({ 0, 0, NULL }),
    _compression({ 0, 0, NULL }),
    _type(other._type) {

    initialize();

    kputsn(other._path.s, other._path.l, &_path);
    kputsn(other._name.s, other._name.l, &_name);
    kputsn(other._directory.s, other._directory.l, &_directory);
    kputsn(other._extension.s, other._extension.l, &_extension);
    kputsn(other._compression.s, other._compression.l, &_compression);
};
URL::URL(const string& path, const IoDirection& direction) : 
    _path({ 0, 0, NULL }),
    _name({ 0, 0, NULL }),
    _directory({ 0, 0, NULL }),
    _extension({ 0, 0, NULL }),
    _compression({ 0, 0, NULL }),
    _type(FormatType::UNKNOWN) {

    initialize();

    parse(path.c_str(), path.size(), direction);
};
URL::URL(const char* path, const size_t size, const IoDirection& direction) : 
    _path({ 0, 0, NULL }),
    _name({ 0, 0, NULL }),
    _directory({ 0, 0, NULL }),
    _extension({ 0, 0, NULL }),
    _compression({ 0, 0, NULL }),
    _type(FormatType::UNKNOWN) {

    initialize();

    parse(path, size, direction);
};
void URL::parse(const char* path, const size_t size, const IoDirection& direction) {
    clear();
    if(size > 0) {
        kputsn(path, size, &_path);

        // first resolve any environment variables in the path
        expand(&_path);

        // standard stream handling 
        // first allow - as an alias for stdin and stdout according to IO direction
        // otherwise check for non canonical paths to standard streams and replace with
        // the canonical one

        if(!strcmp(_path.s, "-")) {
            switch(direction) {
                case IoDirection::IN: {
                    ks_clear(_path);
                    kputs(CANONICAL_STDIN_PATH, &_path);
                    break;
                };
                case IoDirection::OUT: {
                    ks_clear(_path);
                    kputs(CANONICAL_STDOUT_PATH, &_path);
                    break;
                };
            }
        } else {
            if (!strcmp(_path.s, "/dev/stdin") || 
                !strcmp(_path.s, "/dev/fd/0") || 
                !strcmp(_path.s, "/proc/self/fd/0")) {
                ks_clear(_path);
                kputs(CANONICAL_STDIN_PATH, &_path);

            } else if(!strcmp(_path.s, "/dev/stdout") ||
                !strcmp(_path.s, "/dev/fd/1") ||
                !strcmp(_path.s, "/proc/self/fd/1")) {
                ks_clear(_path);
                kputs(CANONICAL_STDOUT_PATH, &_path);

            } else if(!strcmp(_path.s, "/dev/stderr") ||
                !strcmp(_path.s, "/dev/fd/2") ||
                !strcmp(_path.s, "/proc/self/fd/2")) {
                ks_clear(_path);
                kputs(CANONICAL_STDERR_PATH, &_path);
            }
        }

        char* position = NULL;

        // split the path into dirname and basename
        position = strrchr(_path.s, PATH_SEPARATOR);
        if(position != NULL) {
            // the path separator is present
            if(position < _path.s + _path.l) {
                kputs(position + 1, &_name);
            }
            if(position > _path.s) {
                // the last path separator is not the first character
                kputsn(_path.s, position - _path.s, &_directory);
            } else {
                // the last path separator is first character
                // so this is a file in the root directory
                kputc(PATH_SEPARATOR, &_directory);
            }
        } else {
            // relative file path in current directory
            kputsn(_path.s, _path.l, &_name);
        }

        if(!is_standard_stream()) {
            // split extension
            position = strrchr(_name.s, EXTENSION_SEPARATOR);
            if(position != NULL) {
                if(position < _name.s + _name.l) {
                    kputs(position + 1, &_extension);
                }
                _name.l = position - _name.s;
                ks_terminate(_name);

                // if there is a second extension 
                // and the first is a compression marker
                position = strrchr(_name.s, EXTENSION_SEPARATOR);
                if  (position != NULL && (
                    !strcmp(_extension.s, "gz") || 
                    !strcmp(_extension.s, "bz2") || 
                    !strcmp(_extension.s, "xz"))) {

                    kputsn(_extension.s, _extension.l, &_compression);
                    ks_clear(_extension);
                    if(position < _name.s + _name.l) {
                        kputs(position + 1, &_extension);
                    }
                    _name.l = position - _name.s;
                    ks_terminate(_name);
                }
            }

            if(_extension.l > 0 && _type == FormatType::UNKNOWN) {
                _extension.s >> _type;
            }
        }
        refresh();
    }
};
void URL::set_name(const char* name, const size_t size) {
    ks_clear(_name);
    if(name != NULL) {
        kputsn(name, size, &_name);
    }
    refresh();
};
void URL::set_directory(const char* directory, const size_t size) {
    ks_clear(_directory);
    if(directory != NULL) {
        kputsn(directory, size, &_directory);
        expand(&_directory);
    }
    refresh();
};
void URL::set_compression(const char* compression, const size_t size) {
    ks_clear(_compression);
    if(compression != NULL) {
        kputsn(compression, size, &_compression);
    }
    refresh();
};
void URL::set_type(const char* type) {
    type >> _type;
    if(!is_standard_stream()) {
        if(_type != FormatType::UNKNOWN) {
            decode_extension(_type);
        }
        refresh();
    }
};
void URL::set_type(const FormatType type) {
    _type = type;
    if(!is_standard_stream()) {
        if(_type != FormatType::UNKNOWN) {
            decode_extension(_type);
        }
        refresh();
    }
};
void URL::relocate(const URL& base) {
    if(base._directory.l > 0 && !is_absolute()) {
        kstring_t joined = { 0, 0, NULL };
        kputsn(base._directory.s, base._directory.l, &joined);
        if(_directory.l > 0) {
            if(joined.s[joined.l - 1] != PATH_SEPARATOR) {
                kputc(PATH_SEPARATOR, &joined);
            }
            kputsn(_directory.s, _directory.l, &joined);
        }
        kputsn(joined.s, joined.l, &_directory);
        ks_free(joined);
        refresh();
    }
};
FormatKind URL::kind() const {
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
bool URL::is_readable() const {
    if(is_stdin()) {
        return true;
    } else if(is_stdout() || is_stderr() || is_null()) {
        return false;
    } else {
        return access(_path.s, R_OK) != -1;
    }
};
bool URL::is_writable() const {
    if(is_stdin()) {
        return false;
    } else if(is_stdout() || is_stderr() || is_null()) {
        return true;
    } else {
        return access(_path.s, F_OK) != -1 ? access(_path.s, W_OK) != -1 : access(_directory.s, W_OK) != -1;
    }
};
const char* const URL::c_str() const {
    return _path.s;
};
bool URL::operator==(const URL &other) const {
    return !strcmp(other._path.s, _path.s);
};
URL& URL::operator=(const URL& other) {
    if (this != &other) {
        clear();
        kputsn(other._path.s, other._path.l, &_path);
        kputsn(other._name.s, other._name.l, &_name);
        kputsn(other._directory.s, other._directory.l, &_directory);
        kputsn(other._extension.s, other._extension.l, &_extension);
        kputsn(other._compression.s, other._compression.l, &_compression);
        _type = other._type;
    }
    return *this;
};
URL::operator string() const { 
    return string(_path.s, _path.l);
};
void URL::refresh() {
    ks_clear(_path);
    if(_directory.l > 0) {
        kputsn(_directory.s, _directory.l, &_path);
    }
    if(_name.l > 0) {
        if(_path.l > 0 && _path.s[_path.l - 1] != PATH_SEPARATOR) {
            kputc(PATH_SEPARATOR, &_path);
        }
        kputsn(_name.s, _name.l, &_path);
    }
    if(!is_standard_stream()) {
        if(_extension.l > 0) {
            kputc(EXTENSION_SEPARATOR, &_path);
            kputsn(_extension.s, _extension.l, &_path);
        }
        if(_compression.l > 0) {
            switch(_type) {
                case FormatType::SAM:
                case FormatType::FASTQ:
                    kputc(EXTENSION_SEPARATOR, &_path);
                    kputsn(_compression.s, _compression.l, &_path);
                    break;
                default:
                    break;
            }
        }
    }
};
void URL::decode_extension(const FormatType& type) {
    ks_clear(_extension);
    switch(type) {
        case FormatType::FASTQ:
            kputs("fastq", &_extension);
            break;
        case FormatType::SAM:
            kputs("sam", &_extension);
            break;
        case FormatType::BAM:
            kputs("bam", &_extension);
            break;
        case FormatType::BAI:
            kputs("bai", &_extension);
            break;
        case FormatType::CRAM:
            kputs("cram", &_extension);
            break;
        case FormatType::CRAI:
            kputs("crai", &_extension);
            break;
        case FormatType::VCF:
            kputs("vcf", &_extension);
            break;
        case FormatType::BCF:
            kputs("bcf", &_extension);
            break;
        case FormatType::CSI:
            kputs("csi", &_extension);
            break;
        case FormatType::GZI:
            kputs("gzi", &_extension);
            break;
        case FormatType::TBI:
            kputs("tbi", &_extension);
            break;
        case FormatType::BED:
            kputs("bed", &_extension);
            break;
        default:
            break;
    }
};
void URL::encode(Document& document, Value& value) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    if(!is_standard_stream()) {
        value.SetString(_path.s, _path.l, allocator);
    } else {
        value.SetObject();
        Value v;

        v.SetString(_path.s, _path.l, allocator);
        value.AddMember("path", v, allocator);

        if(_compression.l > 0) {
            v.SetString(_compression.s, _compression.l, allocator);
            value.AddMember("compression", v, allocator);
        }

        if(_type != FormatType::UNKNOWN) {
            string type;
            type << _type;
            v.SetString(type.c_str(), type.size(), allocator);
            value.AddMember("type", v, allocator);
        }
    }
};
void URL::expand(kstring_t* path) {
    if(path->l > 0) {
        kstring_t resolved = { 0, 0, NULL };
        kstring_t variable = { 0, 0, NULL };
        char* value = NULL;
        size_t position = 0;
        while(position < path->l) {
            const char& c = path->s[position];
            switch(c) {
                case '~':
                    if(resolved.l == 0) {
                        value = getenv("HOME");
                        if(value != NULL) {
                            kputsn_(value, strlen(value), &resolved);
                        } else {
                            kputc_(c, &resolved);
                        }
                    } else {
                        kputc_(c, &resolved);
                    }
                    break;
                case '$':
                    if(variable.l == 0) {
                        kputc_(c, &variable);
                    } else {
                        kputc_(c, &resolved);
                    }
                    break;
                case '{':
                    if(variable.l == 1) {
                        kputc_(c, &variable);
                    } else {
                        kputc_(c, &resolved);
                    }
                    break;
                case '}':
                    if(variable.l == 0) {
                        kputc_(c, &resolved);
                    } else if(variable.l == 1) {
                        kputc_(*variable.s, &resolved);
                        kputc_(c, &resolved);
                        ks_clear(variable);
                    } else {
                        ks_terminate(variable);
                        value = getenv(variable.s + 2);
                        if(value != NULL) {
                            kputsn_(value, strlen(value), &resolved);
                        }
                        ks_clear(variable);
                    }
                    break;
                case '\0':
                    throw InternalError("unexpected string terminal encountered in " + *this);
                    break;
                default:
                    if(variable.l == 0) {
                        kputc_(c, &resolved);
                    } else if(variable.l == 1) {
                        kputc_(*variable.s, &resolved);
                        kputc_(c, &resolved);
                        ks_clear(variable);
                    } else {
                        kputc_(c, &variable);
                    }
            }
            position++;
        }
        if(variable.l == 0) {
            ks_terminate(resolved);
        } else if(variable.l == 1) {
            kputc(*variable.s, &resolved);
            ks_clear(variable);
        } else {
            string message("unterminated environment variable in path ");
            message += path->s;
            message += " at position ";
            message += to_string(position - variable.l);
            throw ConfigurationError(message);
        }
        ks_clear(*path);
        kputsn(resolved.s, resolved.l, path);
        ks_free(resolved);
        ks_free(variable);
    }
};
string operator+(const string& lhs, const URL& rhs) {
    string value(lhs);
    value.append(rhs._path.s, rhs._path.l);
    return value;
};
string operator+(const URL& lhs, const string& rhs) {
    string value(lhs._path.s, lhs._path.l);
    value += rhs;
    return value;
};
bool operator<(const URL& lhs, const URL& rhs) {
    return lhs._path.s < rhs._path.s;
};
ostream& operator<<(ostream& o, const URL& url) {
    o << url._path.s;
    return o;
};
*/

/*  Token
*/
Token::Token(
    const size_t& index,
    const size_t& input_segment_index,
    const int32_t& start,
    const int32_t& end,
    const bool& end_terminated) :
    index(index),
    input_segment_index(input_segment_index),
    start(start),
    end(end),
    end_terminated(end_terminated) {
};
Token::Token(const Token& other) :
    index(other.index),
    input_segment_index(other.input_segment_index),
    start(other.start),
    end(other.end),
    end_terminated(other.end_terminated) {
};
string Token::description() const {
    string o("cycles ");
    o.append(to_string(start));
    o.append(" to ");
    o.append(end_terminated ? to_string(end) : "end");
    o.append(" of input segment ");
    o.append(to_string(input_segment_index));
    return o;
};
Token::operator string() const {
    string o;
    o.append(to_string(input_segment_index));
    o.push_back(':');
    if(start) o.append(to_string(start));
    o.push_back(':');
    if(end_terminated) o.append(to_string(end));
    return o;
};
ostream& operator<<(ostream& o, const Token& token) {
    o << token.input_segment_index;
    o << ":";
    if(token.start) o << token.start;
    o << ":";
    if(token.end_terminated) o << token.end;
    return o;
};

/*  Transform
*/
Transform::Transform (
    const size_t& index,
    const Token& token,
    const size_t& output_segment_index,
    const LeftTokenOperator& left) :

    index(index),
    output_segment_index(output_segment_index),
    token(token),
    left(left) {
};
Transform::Transform(const Transform& other) :
    index(other.index),
    output_segment_index(other.output_segment_index),
    token(other.token),
    left(other.left) {
};
string Transform::description() const {
    string o("Append ");
    switch (left) {
        case LeftTokenOperator::NONE:
            o.append("token ");
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o.append("reverse complemented token ");
            break;
    }
    o.append(to_string(token.index));
    o.append(" of input segment ");
    o.append(to_string(token.input_segment_index));
    o.append(" to output segment ");
    o.append(to_string(output_segment_index));
    return o;
};
Transform::operator string() const {
    string o;
    switch (left) {
        case LeftTokenOperator::NONE:
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o.push_back('~');
            break;
    }
    o.append(to_string(token.index));
    return o;
};
ostream& operator<<(ostream& o, const Transform& transform) {
    switch (transform.left) {
        case LeftTokenOperator::NONE:
            break;
        case LeftTokenOperator::REVERSE_COMPLEMENT:
            o << '~';
            break;
    }
    o << transform.token.index;
    return o;
};

/*  Sequence
*/
Sequence::Sequence() :
    code(NULL),
    quality(NULL),
    capacity(INITIAL_SEQUENCE_CAPACITY),
    length(0) {
    code = (uint8_t*)malloc(capacity);
    quality = (uint8_t*)malloc(capacity);
    terminate();
};
Sequence::Sequence(const Sequence& other) :
    capacity(other.capacity),
    length(other.length) {
    code = (uint8_t*)malloc(capacity);
    quality = (uint8_t*)malloc(capacity);
    memcpy(code, other.code, length);
    memcpy(quality, other.quality, length);
    terminate();
};
Sequence& Sequence::operator=(const Sequence& other) {
    if(&other == this) {
        return *this;
    } else {
        fill(other.code, other.quality, other.length);
    }
    return *this;
};
Sequence::~Sequence() {
    free(code);
    free(quality);
    code = NULL;
    quality = NULL;
};
void Sequence::mask(const uint8_t& threshold) {
    if (threshold > 0) {
        for (size_t i = 0; i < length; i++) {
            if (quality[i] < threshold) {
                code[i] = ANY_NUCLEOTIDE;
                // TODO do we also set the quality value to 2?
            }
        }
    }
};
void Sequence::expected_error(float& error) const {
    double value = 0.0;
    for (uint8_t* q = quality; *q; ++q) {
        value += quality_to_probability(uint8_t(*q));
    }
    error = float(value);
};
void Sequence::fill(const uint8_t* code, const uint8_t* quality, const size_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }
        memcpy(this->code, code, size);
        memcpy(this->quality, quality, size);
    }
    length = size;
    terminate();
};
void Sequence::fill(const char* code, const size_t& size) {
    if(size > 0) {
        if(size >= capacity) {
            capacity = size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }

        for(size_t i = 0; i < size; i++) {
            *(this->code + i) = AsciiToAmbiguousBam[size_t(code[i])];
            *(this->quality + i) = MAX_VALID_PHRED_VALUE;
        }
    }
    length = size;
    terminate();
};
void Sequence::append(const uint8_t* code, const uint8_t* quality, const size_t& size) {
    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length+ size + 1;
            kroundup32(capacity);
            this->code = (uint8_t*)realloc(this->code, capacity);
            this->quality = (uint8_t*)realloc(this->quality, capacity);
        }
        memcpy(this->code + length, code, size);
        memcpy(this->quality + length, quality, size);
        length += size;
        terminate();
    }
};
void Sequence::append(const Sequence& other, const size_t& start, const size_t& size) {
    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length + size + 1;
            kroundup32(capacity);
            code = (uint8_t*)realloc(code, capacity);
            quality = (uint8_t*)realloc(quality, capacity);
        }
        memcpy(code + length, other.code + start, size);
        memcpy(quality + length, other.quality + start, size);
        length += size;
        terminate();
    }
};
size_t Sequence::append(const Sequence& other, const Transform& transform) {
    const size_t start = transform.token.decode_start(other.length);
    const size_t end = transform.token.decode_end(other.length);
    const size_t size = end - start;

    if(size > 0) {
        if(length + size >= capacity) {
            capacity = length + size + 1;
            kroundup32(capacity);
            code = (uint8_t*)realloc(code, capacity);
            quality = (uint8_t*)realloc(quality, capacity);
        }
        switch (transform.left) {

            case LeftTokenOperator::NONE: {
                memcpy(code + length, other.code + start, size);
                memcpy(quality + length, other.quality + start, size);
                break;
            };

            case LeftTokenOperator::REVERSE_COMPLEMENT: {
                for(size_t i = 0; i < size; i++) {
                    code[length + i] = BamToReverseComplementBam[other.code[end - i - 1]];
                    quality[length + i] = other.quality[end - i - 1];
                }
                break;
            };
        }
        length += size;
        terminate();
    }
    return size;
};
ostream& operator<<(ostream& o, const Sequence& sequence) {
    if(sequence.length > 0) {
        string word;
        sequence.encode_iupac_ambiguity(word);
        word.push_back(LINE_BREAK);
        sequence.encode_phred_quality(word, SAM_PHRED_DECODING_OFFSET);
        word.push_back(LINE_BREAK);
        o << word;
    }
    return o;
};
bool operator<(const Sequence& left, const Sequence& right) {
    size_t position = 0;
    while(position < left.length && position < right.length) {
        if(left.code[position] == right.code[position]) {
            position++;
        } else {
            if(left.code[position] < right.code[position]) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
};
bool operator>(const Sequence& left, const Sequence& right) {
    size_t position = 0;
    while(position < left.length && position < right.length) {
        if(left.code[position] == right.code[position]) {
            position++;
        } else {
            if(left.code[position] > right.code[position]) {
                return true;
            } else {
                return false;
            }
        }
    }
    return false;
};

/*  Barcode
*/
Barcode::Barcode() :
    length(0),
    threshold(0) {
};
Barcode::Barcode(const size_t& width) : 
    length(0),
    tolerance(width),
    fragments(width) {
};
Barcode::Barcode(const Barcode& other) :
    length(other.length),
    threshold(other.threshold),
    tolerance(other.tolerance),
    fragments(other.fragments) {
};
Barcode& Barcode::operator=(const Barcode& other) {
    if(&other == this) {
        return *this;
    } else {
        length = other.length;
        fragments = other.fragments;
        tolerance = other.tolerance;
        threshold = other.threshold;
    }
    return *this;
};
Barcode::operator string() const {
    /* NOTICE barcode is converted to the BAM encoding string, not iupac */
    string key;
    for(const auto& sequence : fragments) {
        for(size_t i = 0; i < sequence.length; i++) {
            key.push_back(sequence.code[i]);
        }
    }
    return key;
};
void Barcode::set_tolerance(const vector<uint8_t>& tolerance) {
    this->tolerance = tolerance;
};
void Barcode::set_threshold(const uint8_t& threshold) {
    this->threshold = threshold;
};
void Barcode::fill(const size_t& position, const char* code, const size_t& size) {
    resize(position + 1);
    length += size;
    fragments[position].fill(code, size);
};
void Barcode::append(const size_t& position, const Sequence& sequence, const Transform& transform) {
    length += fragments[transform.output_segment_index].append(sequence, transform);
};
string Barcode::iupac_ambiguity(const size_t position) const {
    return fragments[position].iupac_ambiguity();
};
string Barcode::iupac_ambiguity() const {
    string result;
    for(const auto& sequence : fragments) {
        for(size_t i = 0; i < sequence.length; i++) {
            result.push_back(BamToAmbiguousAscii[uint8_t(sequence.code[i])]);
        }
    }
    return result;
};
void Barcode::encode_configuration(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Document::AllocatorType& allocator = document.GetAllocator();
        Value code;
        Value collection;
        collection.SetArray();
        for(auto& fragment : fragments) {
            fragment.encode_iupac_ambiguity(code);
            collection.PushBack(code, allocator);
        }
        node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), collection, allocator);
    }
};
void Barcode::encode_report(Document& document, Value& node, const string& key) const {
    if(!empty()) {
        Document::AllocatorType& allocator = document.GetAllocator();
        Value code;
        Value barcode;
        Value collection;
        collection.SetArray();

        for(auto& fragment : fragments) {
            barcode.SetObject();
            fragment.encode_iupac_ambiguity(code);
            barcode.AddMember("barcode sequence", code, allocator);
            collection.PushBack(barcode, allocator);
        }
        node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), collection, allocator);
    }
};
ostream& operator<<(ostream& o, const Barcode& barcode) {
    o << barcode.iupac_ambiguity();
    return o;
};

/* @HD The header line
*/
HeadHDAtom::HeadHDAtom() : 
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }) {
};
HeadHDAtom::HeadHDAtom(const HeadHDAtom& other) :
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }){
    if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
    if(other.SO.l > 0) kputsn(other.SO.s, other.SO.l, &SO);
    if(other.GO.l > 0) kputsn(other.GO.s, other.GO.l, &GO);
};
HeadHDAtom::~HeadHDAtom() {
    ks_free(VN);
    ks_free(SO);
    ks_free(GO);
};
HeadHDAtom& HeadHDAtom::operator=(const HeadHDAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(VN);
        ks_clear(SO);
        ks_clear(GO);
        if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
        if(other.SO.l > 0) kputsn(other.SO.s, other.SO.l, &SO);
        if(other.GO.l > 0) kputsn(other.GO.s, other.GO.l, &GO);
    }
    return *this;
};
void HeadHDAtom::encode(kstring_t* buffer) const {
    kputsn_("@HD", 3, buffer);
    if(VN.l > 0) {
        kputsn_("\tVN:", 4, buffer);
        kputsn_(VN.s, VN.l, buffer);
    }
    if(SO.l > 0) {
        kputsn_("\tSO:", 4, buffer);
        kputsn_(SO.s, SO.l, buffer);
    }
    if(GO.l > 0) {
        kputsn_("\tGO:", 4, buffer);
        kputsn_(GO.s, GO.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadHDAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::VN): {
                position = copy_until_tag_end(position, end, &VN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SO): {
                position = copy_until_tag_end(position, end, &SO);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::GO): {
                position = copy_until_tag_end(position, end, &GO);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadHDAtom::set_alignment_sort_order(const HtsSortOrder& order) {
    SO << order;
};
void HeadHDAtom::set_alignment_grouping(const HtsGrouping& grouping) {
    GO << grouping;
};
void HeadHDAtom::set_version(const htsFormat* format) {
    ks_clear(VN);
    if(format != NULL) {
        kputw(format->version.major, &VN);
        kputc('.', &VN);
        kputw(format->version.minor, &VN);
    }
};
ostream& operator<<(ostream& o, const HeadHDAtom& hd) {
    if(hd.VN.l > 0) o << "VN : " << hd.VN.s << endl;
    if(hd.SO.l > 0) o << "SO : " << hd.SO.s << endl;
    if(hd.GO.l > 0) o << "GO : " << hd.GO.s << endl;
    return o;
};

/* @SQ reference sequence dictionary
*/
HeadSQAtom::HeadSQAtom() :
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){
    ks_terminate(SN)
};
HeadSQAtom::HeadSQAtom(const HeadSQAtom& other) :
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){
    ks_terminate(SN)
    if(other.SN.l > 0) kputsn(other.SN.s, other.SN.l, &SN);
    if(other.LN > 0)   LN = other.LN;
    if(other.AH.l > 0) kputsn(other.AH.s, other.AH.l, &AH);
    if(other.AS.l > 0) kputsn(other.AS.s, other.AS.l, &AS);
    if(other.M5.l > 0) kputsn(other.M5.s, other.M5.l, &M5);
    if(other.SP.l > 0) kputsn(other.SP.s, other.SP.l, &SP);
    if(other.UR.l > 0) kputsn(other.UR.s, other.UR.l, &UR);
};
HeadSQAtom::~HeadSQAtom() {
    ks_free(SN);
    ks_free(AH);
    ks_free(AS);
    ks_free(M5);
    ks_free(SP);
    ks_free(UR);
};
HeadSQAtom& HeadSQAtom::operator=(const HeadSQAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(SN);
        LN = 0;
        ks_clear(AH);
        ks_clear(AS);
        ks_clear(M5);
        ks_clear(SP);
        ks_clear(UR);
        if(other.SN.l > 0) kputsn(other.SN.s, other.SN.l, &SN);
        if(other.LN > 0)   LN = other.LN;
        if(other.AH.l > 0) kputsn(other.AH.s, other.AH.l, &AH);
        if(other.AS.l > 0) kputsn(other.AS.s, other.AS.l, &AS);
        if(other.M5.l > 0) kputsn(other.M5.s, other.M5.l, &M5);
        if(other.SP.l > 0) kputsn(other.SP.s, other.SP.l, &SP);
        if(other.UR.l > 0) kputsn(other.UR.s, other.UR.l, &UR);
    }
    return *this;
};
HeadSQAtom::operator string() const {
    return string(SN.s, SN.l);
};
void HeadSQAtom::encode(kstring_t* buffer) const {
    kputsn_("@SQ", 3, buffer);
    if(SN.l > 0) {
        kputsn_("\tSN:", 4, buffer);
        kputsn_(SN.s, SN.l, buffer);
    }
    if(LN > 0) {
        kputsn_("\tLN:", 4, buffer);
        kputw(LN, buffer);
    }
    if(AH.l > 0) {
        kputsn_("\tAH:", 4, buffer);
        kputsn_(AH.s, AH.l, buffer);
    }
    if(AS.l > 0) {
        kputsn_("\tAS:", 4, buffer);
        kputsn_(AS.s, AS.l, buffer);
    }
    if(M5.l > 0) {
        kputsn_("\tM5:", 4, buffer);
        kputsn_(M5.s, M5.l, buffer);
    }
    if(SP.l > 0) {
        kputsn_("\tSP:", 4, buffer);
        kputsn_(SP.s, SP.l, buffer);
    }
    if(UR.l > 0) {
        kputsn_("\tUR:", 4, buffer);
        kputsn_(UR.s, UR.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadSQAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::SN): {
                position = copy_until_tag_end(position, end, &SN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::LN): {
                LN = strtol(position, &position, 10);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::AH): {
                position = copy_until_tag_end(position, end, &AH);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::AS): {
                position = copy_until_tag_end(position, end, &AS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::M5): {
                position = copy_until_tag_end(position, end, &M5);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SP): {
                position = copy_until_tag_end(position, end, &SP);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::UR): {
                position = copy_until_tag_end(position, end, &UR);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
ostream& operator<<(ostream& o, const HeadSQAtom& sq) {
    if(sq.SN.l > 0) o << "SN : " << sq.SN.s << endl;
    if(sq.LN   > 0) o << "LN : " << sq.LN   << endl;
    if(sq.AH.l > 0) o << "AH : " << sq.AH.s << endl;
    if(sq.AS.l > 0) o << "AS : " << sq.AS.s << endl;
    if(sq.M5.l > 0) o << "M5 : " << sq.M5.s << endl;
    if(sq.SP.l > 0) o << "SP : " << sq.SP.s << endl;
    if(sq.UR.l > 0) o << "UR : " << sq.UR.s << endl;
    return o;
};

/* @PG program
*/
HeadPGAtom::HeadPGAtom() :
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){
    ks_terminate(ID)
};
HeadPGAtom::HeadPGAtom(const HeadPGAtom& other) :
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){
    ks_terminate(ID)
    if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
    if(other.PN.l > 0) kputsn(other.PN.s, other.PN.l, &PN);
    if(other.CL.l > 0) kputsn(other.CL.s, other.CL.l, &CL);
    if(other.PP.l > 0) kputsn(other.PP.s, other.PP.l, &PP);
    if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
    if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
};
HeadPGAtom::~HeadPGAtom() {
    ks_free(ID);
    ks_free(PN);
    ks_free(CL);
    ks_free(PP);
    ks_free(DS);
    ks_free(VN);
};
HeadPGAtom& HeadPGAtom::operator=(const HeadPGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(ID);
        ks_clear(PN);
        ks_clear(CL);
        ks_clear(PP);
        ks_clear(DS);
        ks_clear(VN);
        if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
        if(other.PN.l > 0) kputsn(other.PN.s, other.PN.l, &PN);
        if(other.CL.l > 0) kputsn(other.CL.s, other.CL.l, &CL);
        if(other.PP.l > 0) kputsn(other.PP.s, other.PP.l, &PP);
        if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(other.VN.l > 0) kputsn(other.VN.s, other.VN.l, &VN);
    }
    return *this;
};
HeadPGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadPGAtom::encode(kstring_t* buffer) const {
    kputsn_("@PG", 3, buffer);
    if(ID.l > 0) {
        kputsn_("\tID:", 4, buffer);
        kputsn_(ID.s, ID.l, buffer);
    }
    if(PN.l > 0) {
        kputsn_("\tPN:", 4, buffer);
        kputsn_(PN.s, PN.l, buffer);
    }
    if(CL.l > 0) {
        kputsn_("\tCL:", 4, buffer);
        kputsn_(CL.s, CL.l, buffer);
    }
    if(PP.l > 0) {
        kputsn_("\tPP:", 4, buffer);
        kputsn_(PP.s, PP.l, buffer);
    }
    if(DS.l > 0) {
        kputsn_("\tDS:", 4, buffer);
        kputsn_(DS.s, DS.l, buffer);
    }
    if(VN.l > 0) {
        kputsn_("\tVN:", 4, buffer);
        kputsn_(VN.s, VN.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadPGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::ID): {
                position = copy_until_tag_end(position, end, &ID);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PN): {
                position = copy_until_tag_end(position, end, &PN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::CL): {
                position = copy_until_tag_end(position, end, &CL);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PP): {
                position = copy_until_tag_end(position, end, &PP);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DS): {
                position = copy_until_tag_end(position, end, &DS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::VN): {
                position = copy_until_tag_end(position, end, &VN);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
ostream& operator<<(ostream& o, const HeadPGAtom& pg) {
    if(pg.ID.l > 0) o << "ID : " << pg.ID.s << endl;
    if(pg.PN.l > 0) o << "PN : " << pg.PN.s << endl;
    if(pg.CL.l > 0) o << "CL : " << pg.CL.s << endl;
    if(pg.PP.l > 0) o << "PP : " << pg.PP.s << endl;
    if(pg.DS.l > 0) o << "DS : " << pg.DS.s << endl;
    if(pg.VN.l > 0) o << "VN : " << pg.VN.s << endl;
    return o;
};

/* @RG Read Group
*/
HeadRGAtom::HeadRGAtom() :
    ID({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    SM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }){
    ks_terminate(ID)
};
HeadRGAtom::HeadRGAtom(const HeadRGAtom& other) :
    ID({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    SM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }){
    ks_terminate(ID)
    if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
    if(other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
    if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
    if(other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
    if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
    if(other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
    if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
    if(other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
    if(other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
    if(other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
    if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
    if(other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
    if(other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
};
HeadRGAtom::~HeadRGAtom() {
    ks_free(ID);
    ks_free(PI);
    ks_free(LB);
    ks_free(SM);
    ks_free(PU);
    ks_free(CN);
    ks_free(DS);
    ks_free(DT);
    ks_free(PL);
    ks_free(PM);
    ks_free(PG);
    ks_free(FO);
    ks_free(KS);
};
HeadRGAtom& HeadRGAtom::operator=(const HeadRGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(ID);
        ks_clear(PI);
        ks_clear(LB);
        ks_clear(SM);
        ks_clear(PU);
        ks_clear(CN);
        ks_clear(DS);
        ks_clear(DT);
        ks_clear(PL);
        ks_clear(PM);
        ks_clear(PG);
        ks_clear(FO);
        ks_clear(KS);
        if(other.ID.l > 0) kputsn(other.ID.s, other.ID.l, &ID);
        if(other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
        if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
        if(other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
        if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
        if(other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
        if(other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
        if(other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
        if(other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
        if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
        if(other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
        if(other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
    }
    return *this;
};
HeadRGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadRGAtom::encode(kstring_t* buffer) const {
    kputsn_("@RG", 3, buffer);
    if(ID.l > 0) {
        kputsn_("\tID:", 4, buffer);
        kputsn_(ID.s, ID.l, buffer);
    }
    if(PI.l > 0) {
        kputsn_("\tPI:", 4, buffer);
        kputsn_(PI.s, PI.l, buffer);
    }
    if(LB.l > 0) {
        kputsn_("\tLB:", 4, buffer);
        kputsn_(LB.s, LB.l, buffer);
    }
    if(SM.l > 0) {
        kputsn_("\tSM:", 4, buffer);
        kputsn_(SM.s, SM.l, buffer);
    }
    if(PU.l > 0) {
        kputsn_("\tPU:", 4, buffer);
        kputsn_(PU.s, PU.l, buffer);
    }
    if(CN.l > 0) {
        kputsn_("\tCN:", 4, buffer);
        kputsn_(CN.s, CN.l, buffer);
    }
    if(DS.l > 0) {
        kputsn_("\tDS:", 4, buffer);
        kputsn_(DS.s, DS.l, buffer);
    }
    if(DT.l > 0) {
        kputsn_("\tDT:", 4, buffer);
        kputsn_(DT.s, DT.l, buffer);
    }
    if(PL.l > 0) {
        kputsn_("\tPL:", 4, buffer);
        kputsn_(PL.s, PL.l, buffer);
    }
    if(PM.l > 0) {
        kputsn_("\tPM:", 4, buffer);
        kputsn_(PM.s, PM.l, buffer);
    }
    if(PG.l > 0) {
        kputsn_("\tPG:", 4, buffer);
        kputsn_(PG.s, PG.l, buffer);
    }
    if(FO.l > 0) {
        kputsn_("\tFO:", 4, buffer);
        kputsn_(FO.s, FO.l, buffer);
    }
    if(KS.l > 0) {
        kputsn_("\tKS:", 4, buffer);
        kputsn_(KS.s, KS.l, buffer);
    }
    kputc(LINE_BREAK, buffer);
};
char* HeadRGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        position++;
        uint16_t tag = tag_to_code(position);
        position += 3;
        switch (tag) {
            case uint16_t(HtsAuxiliaryCode::ID): {
                position = copy_until_tag_end(position, end, &ID);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::CN): {
                position = copy_until_tag_end(position, end, &CN);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DS): {
                position = copy_until_tag_end(position, end, &DS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::DT): {
                position = copy_until_tag_end(position, end, &DT);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::FO): {
                position = copy_until_tag_end(position, end, &FO);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::KS): {
                position = copy_until_tag_end(position, end, &KS);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::LB): {
                position = copy_until_tag_end(position, end, &LB);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PG): {
                position = copy_until_tag_end(position, end, &PG);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PI): {
                position = copy_until_tag_end(position, end, &PI);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PL): {
                position = copy_until_tag_end(position, end, &PL);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PM): {
                position = copy_until_tag_end(position, end, &PM);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::PU): {
                position = copy_until_tag_end(position, end, &PU);
                break;
            };
            case uint16_t(HtsAuxiliaryCode::SM): {
                position = copy_until_tag_end(position, end, &SM);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadRGAtom::set_platform(const Platform& value) {
    PL << value;
};
void HeadRGAtom::expand(const HeadRGAtom& other) {
    if(&other != this) {
        if(PI.l == 0 && other.PI.l > 0) kputsn(other.PI.s, other.PI.l, &PI);
        if(LB.l == 0 && other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
        if(SM.l == 0 && other.SM.l > 0) kputsn(other.SM.s, other.SM.l, &SM);
        if(PU.l == 0 && other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
        if(CN.l == 0 && other.CN.l > 0) kputsn(other.CN.s, other.CN.l, &CN);
        if(DS.l == 0 && other.DS.l > 0) kputsn(other.DS.s, other.DS.l, &DS);
        if(DT.l == 0 && other.DT.l > 0) kputsn(other.DT.s, other.DT.l, &DT);
        if(PL.l == 0 && other.PL.l > 0) kputsn(other.PL.s, other.PL.l, &PL);
        if(PM.l == 0 && other.PM.l > 0) kputsn(other.PM.s, other.PM.l, &PM);
        if(PG.l == 0 && other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
        if(FO.l == 0 && other.FO.l > 0) kputsn(other.FO.s, other.FO.l, &FO);
        if(KS.l == 0 && other.KS.l > 0) kputsn(other.KS.s, other.KS.l, &KS);
    }
};
void HeadRGAtom::encode(Document& document, Value& node, const string& key) const {
    Document::AllocatorType& allocator = document.GetAllocator();
    Value v;

    v.SetString(ID.s, ID.l, allocator);
    node.AddMember(Value(key.c_str(), key.size(), allocator).Move(), v, allocator);

    if(PI.l > 0) {
        v.SetString(PI.s, PI.l, allocator);
        node.AddMember("PI", v, allocator);
    }
    if(LB.l > 0) {
        v.SetString(LB.s, LB.l, allocator);
        node.AddMember("LB", v, allocator);
    }
    if(SM.l > 0) {
        v.SetString(SM.s, SM.l, allocator);
        node.AddMember("SM", v, allocator);
    }
    if(PU.l > 0) {
        v.SetString(PU.s, PU.l, allocator);
        node.AddMember("PU", v, allocator);
    }
    if(CN.l > 0) {
        v.SetString(CN.s, CN.l, allocator);
        node.AddMember("CN", v, allocator);
    }
    if(DS.l > 0) {
        v.SetString(DS.s, DS.l, allocator);
        node.AddMember("DS", v, allocator);
    }
    if(DT.l > 0) {
        v.SetString(DT.s, DT.l, allocator);
        node.AddMember("DT", v, allocator);
    }
    if(PL.l > 0) {
        v.SetString(PL.s, PL.l, allocator);
        node.AddMember("PL", v, allocator);
    }
    if(PM.l > 0) {
        v.SetString(PM.s, PM.l, allocator);
        node.AddMember("PM", v, allocator);
    }
    if(PG.l > 0) {
        v.SetString(PG.s, PG.l, allocator);
        node.AddMember("PG", v, allocator);
    }
    if(FO.l > 0) {
        v.SetString(FO.s, FO.l, allocator);
        node.AddMember("FO", v, allocator);
    }
    if(KS.l > 0) {
        v.SetString(KS.s, KS.l, allocator);
        node.AddMember("KS", v, allocator);
    }
};
ostream& operator<<(ostream& o, const HeadRGAtom& rg) {
    if(rg.ID.l > 0) o << "ID : " << rg.ID.s << endl;
    if(rg.PI.l > 0) o << "PI : " << rg.PI.s << endl;
    if(rg.LB.l > 0) o << "LB : " << rg.LB.s << endl;
    if(rg.SM.l > 0) o << "SM : " << rg.SM.s << endl;
    if(rg.PU.l > 0) o << "PU : " << rg.PU.s << endl;
    if(rg.CN.l > 0) o << "CN : " << rg.CN.s << endl;
    if(rg.DS.l > 0) o << "DS : " << rg.DS.s << endl;
    if(rg.DT.l > 0) o << "DT : " << rg.DT.s << endl;
    if(rg.PL.l > 0) o << "PL : " << rg.PL.s << endl;
    if(rg.PM.l > 0) o << "PM : " << rg.PM.s << endl;
    if(rg.PG.l > 0) o << "PG : " << rg.PG.s << endl;
    if(rg.FO.l > 0) o << "FO : " << rg.FO.s << endl;
    if(rg.KS.l > 0) o << "KS : " << rg.KS.s << endl;
    return o;
};

/* @CO free text comment
*/
HeadCOAtom::HeadCOAtom() :
    CO({ 0, 0, NULL }){
};
HeadCOAtom::HeadCOAtom(const HeadCOAtom& other) :
    CO({ 0, 0, NULL }){
    if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
};
HeadCOAtom::~HeadCOAtom() {
    ks_free(CO);
};
HeadCOAtom& HeadCOAtom::operator=(const HeadCOAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(CO);
        if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
    }
    return *this;
};
char* HeadCOAtom::decode(char* position, const char* end) {
    if(*position == '\t' && position <= end) {
        position++;
        position = copy_until_linebreak(position, end, &CO);
    }
    return ++position;
};
void HeadCOAtom::encode(kstring_t* buffer) const {
    if(CO.l > 0) {
        kputsn_("@CO:", 4, buffer);
        kputsn_(CO.s, CO.l, buffer);
        kputc(LINE_BREAK, buffer);
    }
};
ostream& operator<<(ostream& o, const HeadCOAtom& co) {
    if(co.CO.l > 0) o << "CO : " << co.CO.s << endl;
    return o;
};

/*  Feed specification
*/
FeedSpecification::FeedSpecification (
    const IoDirection& direction,
    const size_t& index,
    const URL& url,
    const Platform& platform,
    const uint8_t& phred_offset) :

    direction(direction),
    index(index),
    url(url),
    platform(platform),
    capacity(DEFAULT_FEED_CAPACITY),
    resolution(DEFAULT_FEED_RESOLUTION),
    phred_offset(phred_offset),
    hfile(NULL) {
};
void FeedSpecification::set_capacity(const size_t& capacity) {
    if(capacity != this->capacity) {
        size_t aligned = size_t(capacity / resolution) * resolution;
        if(aligned < capacity) {
            aligned += resolution;
        }
        this->capacity = aligned;
    }
};
void FeedSpecification::set_resolution(const size_t& resolution) {
    if(resolution != this->resolution) {
        size_t aligned = size_t(capacity / resolution) * resolution;
        if(aligned < capacity) {
            aligned += resolution;
        }
        this->resolution = resolution;
        this->capacity = aligned;
    }
};
void FeedSpecification::register_rg(const HeadRGAtom& rg) {
    string key(rg);
    auto record = read_group_by_id.find(key);
    if (record == read_group_by_id.end()) {
        read_group_by_id.emplace(make_pair(key, HeadRGAtom(rg)));
    }
};
void FeedSpecification::register_pg(const HeadPGAtom& pg) {
    string key(pg);
    auto record = program_by_id.find(key);
    if (record == program_by_id.end()) {
        program_by_id.emplace(make_pair(key, HeadPGAtom(pg)));
    }
};
void FeedSpecification::describe(ostream& o) const {
    o << "    ";
    switch (direction) {
        case IoDirection::IN:
            o << "Input";
            break;
        case IoDirection::OUT:
            o << "Output";
            break;
    }
    o << " feed No." << index << endl;
    o << "        Type : " << url.type() << endl;
    // if(strlen(url.compression()) > 0) o << "        Compression : " << url.compression() << endl;
    o << "        Resolution : " << resolution << endl;
    o << "        Phred offset : " << to_string(phred_offset) << endl;
    o << "        Platform : " << platform << endl;
    o << "        Buffer capacity : " << capacity << endl;
    o << "        URL : " << url << endl;

    // if(!program_by_id.empty()) {
    //     o << "\tProgram : " << endl;
    //     for(auto& record : program_by_id) {
    //         o << "\t\t" << record.first << endl;
    //     }
    // }
    // if(!read_group_by_id.empty()) {
    //     o << "\tRead Groups : " << endl;
    //     for(auto& record : read_group_by_id) {
    //         o << "\t\t" << record.first << endl;
    //     }
    // }

    o << endl;
};
void FeedSpecification::probe() {
    /*  Probe input file
        
        Here you can potentially use hfile to probe the file
        and verify file format and potentially examine the first read 
    */
    switch(direction) {
        case IoDirection::IN: {
            hfile = hopen(url.c_str(), "r");
            if(url.type() == FormatType::UNKNOWN) {
                size_t buffer_capacity = PEEK_BUFFER_CAPACITY;
                ssize_t buffer_length = 0;
                unsigned char* buffer = (unsigned char*)malloc(buffer_capacity);;

                htsFormat format;
                if(!hts_detect_format(hfile, &format)) {
                    switch (format.format) {
                        case htsExactFormat::sam:
                            url.set_type(FormatType::SAM);
                            break;
                        case htsExactFormat::bam:
                            url.set_type(FormatType::BAM);
                            break;
                        case htsExactFormat::bai:
                            url.set_type(FormatType::BAI);
                            break;
                        case htsExactFormat::cram:
                            url.set_type(FormatType::CRAM);
                            break;
                        case htsExactFormat::crai:
                            url.set_type(FormatType::CRAI);
                            break;
                        case htsExactFormat::vcf:
                            url.set_type(FormatType::VCF);
                            break;
                        case htsExactFormat::bcf:
                            url.set_type(FormatType::BCF);
                            break;
                        case htsExactFormat::csi:
                            url.set_type(FormatType::CSI);
                            break;
                        case htsExactFormat::gzi:
                            url.set_type(FormatType::GZI);
                            break;
                        case htsExactFormat::tbi:
                            url.set_type(FormatType::TBI);
                            break;
                        case htsExactFormat::bed:
                            url.set_type(FormatType::BED);
                            break;
                        default:
                            url.set_type(FormatType::UNKNOWN);
                            break;
                    }
                }

                if(url.type() == FormatType::SAM) {
                    buffer_length = hpeek(hfile, buffer, buffer_capacity);
                    if(buffer_length > 0) {
                        switch (format.compression) {
                            case htsCompression::gzip:
                            case htsCompression::bgzf: {
                                unsigned char* decompressed_buffer = (unsigned char*)malloc(buffer_capacity);;
                                z_stream zstream;
                                zstream.zalloc = NULL;
                                zstream.zfree = NULL;
                                zstream.next_in = buffer;
                                zstream.avail_in = buffer_length;
                                zstream.next_out = decompressed_buffer;
                                zstream.avail_out = buffer_capacity;
                                if(inflateInit2(&zstream, 31) == Z_OK) {
                                    while(zstream.total_out < buffer_capacity) {
                                        if(inflate(&zstream, Z_SYNC_FLUSH) != Z_OK) break;
                                    }
                                    inflateEnd(&zstream);
                                    memcpy(buffer, decompressed_buffer, zstream.total_out);
                                    buffer_length = zstream.total_out;
                                } else {
                                    buffer_length = 0;
                                }
                                free(decompressed_buffer);
                                break;
                            };
                            case htsCompression::no_compression:
                                break;
                            default:
                                throw InternalError("unknown compression");
                                break;
                        }
                    }
                    if(buffer_length > 0) { 
                        size_t state = 0;
                        char* position = (char*)buffer;
                        char* end = position + buffer_length;
                        while(position < end && position != NULL) {
                            if(state == 0) {
                                if(*position == '\n') {
                                    position++;
                                } else {
                                    if(*position == '@') {
                                        state = 1;
                                    } else {
                                        break;
                                    }
                                }
                            } else if(state == 1) {
                                if((*position >= 'A' && *position <= 'Z') || (*position >= 'a' && *position <= 'z')) {
                                    state = 2;
                                } else {
                                    break;
                                }
                            } else if(state == 2) {
                                if(*position == '+' && position < end && *(position + 1) == '\n') {
                                    url.set_type(FormatType::FASTQ);
                                }
                                break;
                            }
                            if((position = strchr(position, '\n')) != NULL) position++;
                        }
                    }
                }
                free(buffer);
            }
            break;
        };
        case IoDirection::OUT: {
            hfile = hopen(url.c_str(), "w");
            break;
        };
    }
};
ostream& operator<<(ostream& o, const FeedSpecification& specification) {
    o << specification.url;
    return o;
};

/*  Channel specification */
ChannelSpecification::ChannelSpecification(size_t index) :
    index(index),
    TC(numeric_limits<size_t>::max()),
    FS({ 0, 0, NULL }),
    CO({ 0, 0, NULL }),
    disable_quality_control(false),
    long_read(false),
    include_filtered(true),
    undetermined(false),
    concentration(0) {
};
ChannelSpecification::~ChannelSpecification() {
    ks_free(FS);
    ks_free(CO);
};
string ChannelSpecification::alias() const {
    string alias("Channel No.");
    alias.append(to_string(index));
    return alias;
};
void ChannelSpecification::describe(ostream& o) const {
    o << "    " << alias() << endl;
    // o << "\tRG index : " << head_read_group->index << endl;
    o << "        RG : " << rg.ID.s << endl;
    if(rg.PU.l > 0) o << "        PU : " << rg.PU.s << endl;
    if(rg.LB.l > 0) o << "        LB : " << rg.LB.s << endl;
    if(rg.SM.l > 0) o << "        SM : " << rg.SM.s << endl;
    if(rg.CN.l > 0) o << "        CN : " << rg.CN.s << endl;
    if(rg.DS.l > 0) o << "        DS : " << rg.DS.s << endl;
    if(rg.DT.l > 0) o << "        DT : " << rg.DT.s << endl;
    if(rg.PL.l > 0) o << "        PL : " << rg.PL.s << endl;
    if(rg.PM.l > 0) o << "        PM : " << rg.PM.s << endl;
    if(rg.PG.l > 0) o << "        PG : " << rg.PG.s << endl;
    if(rg.FO.l > 0) o << "        FO : " << rg.FO.s << endl;
    if(rg.KS.l > 0) o << "        KS : " << rg.KS.s << endl;
    if(rg.PI.l > 0) o << "        PI : " << rg.PI.s << endl;
    if(TC)          o << "        TC : " << TC      << endl;
    if(FS.l > 0)    o << "        FS : " << FS.s    << endl;
    if(CO.l > 0)    o << "        CO : " << CO.s    << endl;
    if(concentration > 0) {
                    o << "        PC : " << setprecision(numeric_limits<double>::digits10 + 1) << concentration << endl;
    }
    if(undetermined) {
        o << "        Undetermined : true" << endl;
    }
    if(!undetermined) {
        o << endl;
        for(size_t i = 0; i < multiplex_barcode.total_fragments(); i++) {
            o << "        Multiplex barcode No." << i << " : " << multiplex_barcode.iupac_ambiguity(i) << endl;
        }
    }
    if (!output_urls.empty()) {
        o << endl;
        for (size_t i = 0; i < output_urls.size(); i++) {
            o << "        Output segment No." << i << " : " << output_urls[i] << endl;
        }
    }
    o << endl;
};
void ChannelSpecification::encode(Document& document, Value& node) const {
    Document::AllocatorType& allocator = document.GetAllocator();

    Value v;
    Value channel;
    channel.SetObject();

    rg.encode(document, channel, "RG");
    if(FS.l > 0) {
        v.SetString(FS.s, FS.l, allocator);
        channel.AddMember("FS", v, allocator);
    }
    if(CO.l > 0) {
        v.SetString(CO.s, CO.l, allocator);
        channel.AddMember("CO", v, allocator);
    }
    if(undetermined) {
        v.SetBool(undetermined);
        channel.AddMember("undetermined", v, allocator);
    } else {
        v.SetDouble(concentration);
        channel.AddMember("concentration", v, allocator);
        multiplex_barcode.encode_configuration(document, channel, "barcode");
    }
    // if(!output_urls.empty()) {
    //     Value collection;
    //     collection.SetArray();
    //     for(auto& url : output_urls) {
    //         url.encode(document, v);
    //         collection.PushBack(v, allocator);
    //     }
    //     channel.AddMember("output", collection, allocator);
    // }
    node.PushBack(channel, allocator);
};
ostream& operator<<(ostream& o, const ChannelSpecification& specification) {
    o << specification.alias();
    return o;
};