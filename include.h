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

#ifndef PHENIQS_INCLUDE_H
#define PHENIQS_INCLUDE_H

/*  Support for accouting for the illumina control number in FASTQ comment field
    #define PHENIQS_ILLUMINA_CONTROL_NUMBER
*/

/*  Support for propegating SAM tags pheniqs does not directly interact with
    #define PHENIQS_EXTENDED_SAM_TAG
*/

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <errno.h>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <mutex>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

using std::cerr;
using std::condition_variable;
using std::cout;
using std::endl;
using std::exception;
using std::fixed;
using std::hash;
using std::ifstream;
using std::int16_t;
using std::int32_t;
using std::int64_t;
using std::int8_t;
using std::invalid_argument;
using std::ios_base;
using std::istreambuf_iterator;
using std::list;
using std::lock_guard;
using std::log10;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::move;
using std::mutex;
using std::numeric_limits;
using std::ostream;
using std::out_of_range;
using std::pair;
using std::recursive_mutex;
using std::set;
using std::setprecision;
using std::setw;
using std::size_t;
using std::string;
using std::thread;
using std::to_string;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;
using std::uint8_t;
using std::unique_lock;
using std::unordered_map;
using std::vector;

#endif /* PHENIQS_INCLUDE_H */
