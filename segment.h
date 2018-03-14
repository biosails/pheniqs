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

#ifndef PHENIQS_SEGMENT_H
#define PHENIQS_SEGMENT_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <htslib/sam.h>
#include <htslib/cram.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/hfile.h>
#include <htslib/thread_pool.h>

#include "error.h"
#include "json.h"
#include "constant.h"
#include "sequence.h"
#include "auxiliary.h"

using std::map;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;

using std::mutex;
using std::recursive_mutex;
using std::condition_variable;
using std::unique_lock;
using std::lock_guard;
using std::thread;

class Segment {
    friend ostream& operator<<(ostream& o, const Segment& segment);
    void operator=(Segment const &) = delete;

    public:
        const size_t index;
        const Platform platform;
        kstring_t name;
        uint16_t flag;
        Sequence sequence;
        Auxiliary auxiliary;
        Segment(const Platform& platform);
        Segment(const size_t& index, const int32_t& FI, const int32_t& TC, const Platform& platform);
        Segment(const Segment& other);
        ~Segment();
        inline int32_t get_segment_index() const {
            if(!auxiliary.FI) {
                if(flag & uint16_t(HtsFlag::PAIRED)) {
                    if(flag & uint16_t(HtsFlag::READ1)) {
                        return 1;
                    } else if(flag & uint16_t(HtsFlag::READ2)) {
                        return 2;
                    } else {
                        throw SequenceError("inconsistent SAM flags");
                    }
                } else {
                    return 1;
                }
            } else {
                return auxiliary.FI;
            }
        };
        inline int32_t get_total_segments() const {
            if(!auxiliary.TC) {
                if(flag & uint16_t(HtsFlag::PAIRED)) {
                    return 2;
                } else {
                    return 1;
                }
            } else {
                return auxiliary.TC;
            }
        };
        inline void set_qcfail(const bool value) {
            if (value) {
                flag |= uint16_t(HtsFlag::QCFAIL);
            } else {
                flag &= ~uint16_t(HtsFlag::QCFAIL);
            }
        };
        inline bool get_qcfail() const {
            return flag & uint16_t(HtsFlag::QCFAIL);
        };
        inline void clear() {
            ks_clear(name);
            set_qcfail(false);
            sequence.clear();
            auxiliary.clear();
        };
};
#endif /* PHENIQS_SEGMENT_H */
