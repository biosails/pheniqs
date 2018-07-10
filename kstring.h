/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    Ported from htslib/kstring.h to be more C++ friendly


    The MIT License

    Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef PHENIQS_KSTRING_H
#define PHENIQS_KSTRING_H

#include "include.h"
#include "error.h"
#include <htslib/kstring.h>

#define tag_to_code(t) static_cast< uint16_t >(*(t)) << 8 | static_cast< uint8_t >(*((t) + 1))
const char LINE_BREAK('\n');

static inline void ks_increase_by_size(kstring_t& s, size_t size) {
    size += s.l;
    if(size > s.m) {
        kroundup_size_t(size);
        char* temp;
        if((temp = static_cast< char* >(realloc(s.s, size))) == NULL) {
            throw OutOfMemoryError();
        }
        s.s = temp;
        s.m = size;
    }
};
static inline void ks_increase_to_size(kstring_t& s, size_t size) {
    if(s.m < size) {
        kroundup_size_t(size);
        char* temp;
        if((temp = static_cast< char* >(realloc(s.s, size))) == NULL) {
            throw OutOfMemoryError();
        }
        s.s = temp;
        s.m = size;
    }
};
static inline bool ks_empty(const kstring_t& s) {
    return s.l == 0;
};
static inline void ks_clear(kstring_t& s) {
    s.l = 0;
    if(s.s != NULL) {
        s.s[0] = '\0';
    }
};
static inline void ks_free(kstring_t& s) {
    if(s.s != NULL) {
        free(s.s);
        s.s = NULL;
    }
    s.l = 0;
    s.m = 0;
};
static inline void ks_terminate(kstring_t& s) {
    ks_increase_by_size(s, 1);
    s.s[s.l] = '\0';
};
static inline void ks_put_character(const char c, kstring_t& s) {
    ks_increase_by_size(s, 2);
    s.s[s.l] = c;
    ++s.l;
    s.s[s.l] = '\0';
};
static inline void ks_put_character_(const char c, kstring_t& s) {
    ks_increase_by_size(s, 1);
    s.s[s.l] = c;
    ++s.l;
};
static inline void ks_put_string(const char* p, size_t l, kstring_t& s) {
    ks_increase_by_size(s, l + 2);
    memcpy(s.s + s.l, p, l);
    s.l += l;
    s.s[s.l] = '\0';
};
static inline void ks_put_string(const kstring_t& from, kstring_t& to) {
    ks_increase_by_size(to, from.l + 2);
    memcpy(to.s + to.l, from.s, from.l);
    to.l += from.l;
    to.s[to.l] = '\0';
};
static inline void ks_put_string(const char* p, kstring_t& s) {
    return ks_put_string(p, strlen(p), s);
};
static inline void ks_put_string(const string& from, kstring_t& to) {
    ks_increase_by_size(to, from.size() + 2);
    memcpy(to.s + to.l, from.c_str(), from.size());
    to.l += from.size();
    to.s[to.l] = '\0';
};
static inline void ks_put_string_(const void* p, size_t l, kstring_t& s) {
    ks_increase_by_size(s, l);
    memcpy(s.s + s.l, p, l);
    s.l += l;
};
static inline void ks_put_string_(const kstring_t& from, kstring_t& to) {
    ks_increase_by_size(to, from.l);
    memcpy(to.s + to.l, from.s, from.l);
    to.l += from.l;
};
static inline void ks_put_uint16(const uint16_t& c, kstring_t& s) {
    char buffer[8];
    uint16_t x;
    int16_t l(0);
    int16_t i(0);

    if(c == 0) {
        ks_put_character('0', s);
    } else {
        for(l = 0, x = c; x > 0; x /= 10) {
            buffer[l] = x % 10 + '0';
            ++l;
        }
        ks_increase_by_size(s, l + 2);
        for(i = l - 1; i >= 0; --i) {
            s.s[s.l] = buffer[i];
            ++s.l;
        }
        s.s[s.l] = '\0';
    }
};
static inline void ks_put_int32(const int32_t& c, kstring_t& s) {
    char buffer[16];
    uint32_t x(c);
    int16_t i(0);
    int16_t l(0);

    if(c < 0) {
        x = -x;
    }
    do {
        buffer[l] = x % 10 + '0';
        ++l;
        x /= 10;
    } while(x > 0);
    if(c < 0) {
        buffer[l] = '-';
        ++l;
    }
    ks_increase_by_size(s, l + 2);
    for(i = l - 1; i >= 0; --i) {
        s.s[s.l] = buffer[i];
        ++s.l;
    }
    s.s[s.l] = '\0';
};
static inline void ks_put_uint32(const uint32_t c, kstring_t& s) {
    char buffer[16];
    uint32_t x;
    int16_t l(0);
    int16_t i(0);

    if(c == 0) {
        ks_put_character('0', s);
    } else {
        for(l = 0, x = c; x > 0; x /= 10) {
            buffer[l] = x % 10 + '0';
            ++l;
        }
        ks_increase_by_size(s, l + 2);
        for(i = l - 1; i >= 0; --i) {
            s.s[s.l] = buffer[i];
            ++s.l;
        }
        s.s[s.l] = '\0';
    }
};
static inline void ks_put_int64(const int64_t c, kstring_t& s) {
    char buffer[32];
    int64_t x(c);
    int i(0);
    int l(0);

    if(c < 0) {
        x = -x;
    }
    do {
        buffer[l] = x % 10 + '0';
        ++l;
        x /= 10;
    } while(x > 0);
    if(c < 0) {
        buffer[l] = '-';
        ++l;
    }
    ks_increase_by_size(s, l + 2);
    for(i = l - 1; i >= 0; --i) {
        s.s[s.l] = buffer[i];
        ++s.l;
    }
    s.s[s.l] = '\0';
};

#endif /* PHENIQS_KSTRING_H */
