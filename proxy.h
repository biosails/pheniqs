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

#ifndef PHENIQS_PROXY_H
#define PHENIQS_PROXY_H

#include "include.h"
#include "url.h"
#include "atom.h"

const ssize_t PEEK_BUFFER_CAPACITY(4096);
const int DEFAULT_FEED_CAPACITY(60);
const int DEFAULT_FEED_RESOLUTION(60);

class FeedProxy {
    friend ostream& operator<<(ostream& o, const FeedProxy& proxy);

    public:
        int32_t index;
        URL url;
        IoDirection direction;
        uint8_t phred_offset;
        hFILE* hfile;
        int capacity;
        int resolution;
        Platform platform;
        unordered_map< string, const HeadPGAtom > program_by_id;
        unordered_map< string, const HeadRGAtom > read_group_by_id;
        FeedProxy(const Value& ontology);
        inline bool is_dev_null() const {
            return url.is_dev_null();
        };
        inline bool is_stdin() const {
            return url.is_stdin();
        };
        inline bool is_stdout() const {
            return url.is_stdout();
        };
        inline bool is_stderr() const {
            return url.is_stderr();
        };
        inline FormatKind kind() const {
            return url.kind();
        };
        void register_rg(const HeadRGAtom& rg);
        void register_pg(const HeadPGAtom& pg);
        void open();
};
bool encode_key_value(const string& key, const FeedProxy& value, Value& container, Document& document);
ostream& operator<<(ostream& o, const FeedProxy& proxy);
template<> FeedProxy decode_value< FeedProxy >(const Value& container);
template<> list< FeedProxy > decode_value_by_key< list< FeedProxy > >(const Value::Ch* key, const Value& container);

#endif /* PHENIQS_PROXY_H */
