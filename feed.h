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

#ifndef PHENIQS_FEED_H
#define PHENIQS_FEED_H

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

#include "constant.h"
#include "error.h"
#include "model.h"

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

KSEQ_INIT(BGZF*, bgzf_read)

#define probability_to_quality(p) uint8_t(-10.0 * log10(p))
#define probability_to_phred(p,o) char(uint8_t(-10.0 * log10(p)) + (o))
#define phred_to_quality(c,o) uint8_t((c) - (o))
#define phred_to_quality_double(c,o) double((c) - (o))
#define quality_to_phred(c,o) char((c) + (o))
#define bam1_seq_seti(s, i, c) ((s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))
#define kseq_clear(x) ks_clear((x)->name); ks_clear((x)->comment); ks_clear((x)->seq); ks_clear((x)->qual);

/*  HTS header
*/
class HtsHeader {
    friend ostream& operator<<(ostream& o, const HtsHeader& head);
    HtsHeader(HtsHeader const &) = delete;
    void operator=(HtsHeader const &) = delete;

    public:
        bam_hdr_t* hdr;
        HeadHDAtom hd;
        vector< HeadCOAtom > comments;
        map< string, const HeadPGAtom > program_by_id;
        map< string, const HeadRGAtom > read_group_by_id;
        HtsHeader();
        ~HtsHeader();
        void assemble();
        void decode(htsFile* hts_file);
        void encode(htsFile* hts_file) const;
        void add_read_group(const HeadRGAtom& read_group);
        void add_program(const HeadPGAtom& program);
        void add_comment(const HeadCOAtom& co);
};
/*  Optional auxiliary fields

    FI  i   The index of segment in the template.
    TC  i   The number of segments in the template.
    RG  Z   Read group. Value matches the header RG-ID tag if @RG is present in the header.
    BC  Z   Multiplex barcode sequence, with any quality scores stored in the QT tag.
    QT  Z   Phred encoded quality of the multiplex barcode sequence in the BC tag.
  * RX  Z   Raw sequence bases of the molecular barcode
  * QX  Z   Raw sequence quality of the molecular barcode
  * BX  Z   Corrected sequence bases of the molecular barcode
  * PX  f   Molecular barcode correction error probability
  * DQ  f   The probability that the demultiplexing decision was incorrect
  * EE  f   The expected number of errors in the segment
    FS  Z   Segment suffix.
    LB  Z   Library. Value to be consistent with the header RG-LB tag if @RG is present.
    PG  Z   Program. Value matches the header PG-ID tag if @PG is present.
    PU  Z   Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
    CO  Z   Free-text comments


    Benchmark mode user space tags

  * XI  i   Illumina control flags

  * XM  Z   MDD barcode
  * YD  i   MDD multiplex distance
  * XL  Z   PAMLD barcode
  * XD  i   PAMLD multiplex distance
  * XP  f   PAMLD conditioned probability
*/
class Auxiliary {
    friend ostream& operator<<(ostream& o, const Auxiliary& auxiliary);
    void operator=(Auxiliary const &) = delete;

    public:
        int32_t FI;
        int32_t TC;
        kstring_t RG;
        kstring_t BC;
        kstring_t QT;
        kstring_t FS;
        kstring_t LB;
        kstring_t PG;
        kstring_t PU;
        kstring_t CO;
        kstring_t RX;
        kstring_t QX;
        kstring_t BX;
        float PX;
        float DQ;
        float EE;
        int32_t XI;
        int32_t YD;
        int32_t XD;
        kstring_t XM;
        kstring_t XL;
        float XP;
        Auxiliary(const int32_t& FI, const int32_t& TC);
        Auxiliary(const Auxiliary& other);
        ~Auxiliary();
        void decode(const bam1_t* bam1);
        void encode(bam1_t* bam1) const;
        inline void set_multiplex_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(&BC);
            barcode.encode_phred_quality(&QT, SAM_PHRED_DECODING_OFFSET);
        };
        inline void set_molecular_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(&RX);
            barcode.encode_phred_quality(&QX, SAM_PHRED_DECODING_OFFSET);
        };
        inline void set_mdd_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(&XM);
        };
        inline void set_pamld_barcode(const Barcode& barcode) {
            barcode.encode_iupac_ambiguity(&XL);
        };
        inline void clear() {
            /* FI and TC don't change */
            ks_clear(RG);
            ks_clear(BC);
            ks_clear(QT);
            ks_clear(FS);
            ks_clear(LB);
            ks_clear(PG);
            ks_clear(PU);
            ks_clear(CO);
            ks_clear(RX);
            ks_clear(QX);
            ks_clear(BX);
            PX = 0;
            DQ = 0;
            EE = 0;
            XI = 0;
            YD = 0;
            XD = 0;
            ks_clear(XM);
            ks_clear(XL);
            XP = 0;
        };
        void imitate(const Auxiliary& other) {
            ks_clear(RG);
            ks_clear(BC);
            ks_clear(QT);
            ks_clear(FS);
            ks_clear(LB);
            ks_clear(PG);
            ks_clear(PU);
            ks_clear(CO);
            ks_clear(RX);
            ks_clear(QX);
            ks_clear(BX);
            if(other.RG.l > 0) kputsn(other.RG.s, other.RG.l, &RG);
            if(other.BC.l > 0) kputsn(other.BC.s, other.BC.l, &BC);
            if(other.QT.l > 0) kputsn(other.QT.s, other.QT.l, &QT);
            if(other.FS.l > 0) kputsn(other.FS.s, other.FS.l, &FS);
            if(other.LB.l > 0) kputsn(other.LB.s, other.LB.l, &LB);
            if(other.PG.l > 0) kputsn(other.PG.s, other.PG.l, &PG);
            if(other.PU.l > 0) kputsn(other.PU.s, other.PU.l, &PU);
            if(other.CO.l > 0) kputsn(other.CO.s, other.CO.l, &CO);
            if(other.RX.l > 0) kputsn(other.RX.s, other.RX.l, &RX);
            if(other.QX.l > 0) kputsn(other.QX.s, other.QX.l, &QX);
            if(other.BX.l > 0) kputsn(other.BX.s, other.BX.l, &BX);
            DQ = other.DQ;
        };
};
/* Read segment
*/
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
/* FASTQ Record
*/
class FastqRecord {
    FastqRecord(FastqRecord const &) = delete;
    void operator=(FastqRecord const &) = delete;

    public:
        kstring_t sequence;
        kstring_t quality;
        kstring_t name;
        kstring_t comment;
        FastqRecord() :
            sequence({ 0, 0, NULL }),
            quality({ 0, 0, NULL }),
            name({ 0, 0, NULL }),
            comment({ 0, 0, NULL }) {
            ks_terminate(sequence);
            ks_terminate(quality);
            ks_terminate(name);
            ks_terminate(comment);
            clear();
        };
        ~FastqRecord() {
            ks_free(sequence);
            ks_free(quality);
            ks_free(name);
            ks_free(comment);
        };
        inline void decode(const kseq_t* kseq, const uint8_t phred_offset) {
            // read a kseq record and populate the FastqRecord
            clear();

            // copy from kseq_t to record
            kputsn(kseq->name.s, kseq->name.l, &name);
            kputsn(kseq->comment.s, kseq->comment.l, &comment);

            // decode sequence
            ks_resize(&sequence, kseq->seq.l + 1);
            for(size_t i = 0; i < kseq->seq.l; i++) {
                sequence.s[i] = AsciiToAmbiguousBam[uint8_t(kseq->seq.s[i])];
                // kputc(AsciiToAmbiguousBam[uint8_t(kseq->seq.s[i])], &sequence);
            }
            sequence.l = kseq->seq.l;
            sequence.s[sequence.l] = '\0';
            // ks_terminate(sequence);

            // decode quality
            ks_resize(&quality, kseq->qual.l + 1);
            for(size_t i = 0; i < kseq->qual.l; i++) {
                quality.s[i] = kseq->qual.s[i] - phred_offset;
                // kputc(kseq->qual.s[i] - phred_offset, &quality);
            }
            sequence.l = kseq->qual.l;
            sequence.s[sequence.l] = '\0';
            // ks_terminate(sequence);
        };
        inline void decode(const Segment& segment) {
            clear();

            // copy from segment to record
            kputsn((char*)segment.sequence.code, segment.sequence.length, &sequence);
            kputsn((char*)segment.sequence.quality, segment.sequence.length, &quality);
            kputsn(segment.name.s, segment.name.l, &name);
            decode_comment(segment);
        };
        inline void encode(Segment& segment) const {

            // write the FastqRecord to Segment
            segment.sequence.fill((const uint8_t*)(sequence.s), (const uint8_t*)(quality.s), sequence.l);
            kputsn(name.s, name.l, &segment.name);
            kputsn(comment.s, comment.l, &segment.auxiliary.CO);
            segment.auxiliary.FI = 0;
            segment.set_qcfail(false);
            segment.auxiliary.XI = 0;
            switch (segment.platform) {
                case Platform::ILLUMINA: {
                    /*  @HWI-ST911:232:HABDFADXX:1:1101:1224:1932 1:N:0:CGATGT
                        name    0:1:2:3:4:5:6
                        0   Instrument ID   HWI-ST911
                        1   Run count       232
                        2   Flowcell ID     HABDFADXX
                        3   Lane number     1
                        4   Tile number     1101
                        5   x coordinate    1224
                        6   y coordinate    1932

                        comment 7:8:9:10
                        7   Segment number  1               uint8_t
                        8   Filtered        N               N|Y
                        9   control number  0               uint16_t
                        10  Barcode         CGATGT          char*
                    */
                    size_t offset = 0;
                    parse_illumina_segment_index(segment, offset);
                    parse_illumina_filtered(segment, offset);
                    parse_illumina_control(segment, offset);
                    parse_illumina_barcode(segment, offset);
                    break;
                }
                case Platform::CAPILLARY:
                case Platform::LS454:
                case Platform::HELICOS:
                case Platform::ONT:
                case Platform::PACBIO:
                case Platform::SOLID:
                case Platform::IONTORRENT:
                    break;
                default:
                    break;
            };
        };
        inline void encode(kstring_t& buffer, uint8_t& phred_offset) const {
            // encode identifier
            kputc('@', &buffer);
            kputsn(name.s, name.l, &buffer);
            kputc(' ', &buffer);
            kputsn(comment.s, comment.l, &buffer);
            kputc(LINE_BREAK, &buffer);

            // encode sequence
            ks_resize(&buffer, buffer.l + sequence.l + 1);
            for (size_t i = 0; i < sequence.l; i++) {
                buffer.s[buffer.l + i] = BamToAmbiguousAscii[uint8_t(sequence.s[i])];
                // kputc(BamToAmbiguousAscii[uint8_t(sequence.s[i])], &buffer);
            }
            buffer.l += sequence.l;
            kputc(LINE_BREAK, &buffer);

            // encode separator
            kputc('+', &buffer);
            kputc(LINE_BREAK, &buffer);

            // encode quality
            ks_resize(&buffer, buffer.l + quality.l + 1);
            for (size_t i = 0; i < quality.l; i++) {
                buffer.s[buffer.l + i] = quality.s[i] + phred_offset;
                // kputc(quality.s[i] + phred_offset, &buffer);
            }
            buffer.l += quality.l;
            kputc(LINE_BREAK, &buffer);
        };

    private:
        inline void clear() {
            ks_clear(sequence);
            ks_clear(quality);
            ks_clear(name);
            ks_clear(comment);
        };
        inline void decode_comment(const Segment& segment) {
            switch (segment.platform) {
                case Platform::CAPILLARY:
                    break;
                case Platform::LS454:
                    break;
                case Platform::ILLUMINA: {
                    kputuw(segment.auxiliary.FI, &comment);
                    kputc(':', &comment);
                    kputc(segment.get_qcfail()? 'Y' : 'N', &comment);
                    kputc(':', &comment);
                    kputuw(segment.auxiliary.XI, &comment);
                    kputc(':', &comment);
                    kputsn(segment.auxiliary.BC.s, segment.auxiliary.BC.l, &comment);
                    break;
                };
                case Platform::SOLID:
                    break;
                case Platform::HELICOS:
                    break;
                case Platform::IONTORRENT:
                    break;
                case Platform::ONT:
                    break;
                case Platform::PACBIO:
                    break;
                default:
                    break;
            };
        };
        inline void parse_illumina_segment_index(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            size_t value = 0;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    // Done parsing this token
                    break;

                } else {
                    if (code >= -1) {
                        if (c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if (code < 0) {
                value = 1;
            }
            segment.auxiliary.FI = value;
            offset = position;
        };
        inline void parse_illumina_filtered(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    // Done parsing this token
                    break;

                } else {
                    if (code == -1) {
                        switch(c) {
                            case 'Y':
                                segment.set_qcfail(true);
                                code = 0;
                                break;
                            case 'N':
                                segment.set_qcfail(false);
                                code = 0;
                                break;
                            default:
                                segment.set_qcfail(false);
                                code = -3;
                                break;
                        }
                    } else {
                        // code equals -2 means there was more than one character to the next separator
                        code = -2;
                    }
                }
            }
            offset = position;
        };
        inline void parse_illumina_control(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            uint16_t value = 0;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    break;

                } else {
                    if (code >= -1) {
                        if (c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if (code < 0) {
                value = 0;
            }
            segment.auxiliary.XI = value;
            offset = position;
        };
        inline void parse_illumina_barcode(Segment& segment, size_t& offset) const {
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ' ') {
                    break;
                }
            }
            size_t length = position - offset;
            if(length > 0) {
                kputsn(comment + offset, length, &(segment.auxiliary.BC));
            }
            offset = position;
        };
};
/* IO feed
*/
class Feed {
    friend class Factory;

    public:
        FeedSpecification const * specification;
        const IoDirection direction;
        const URL url;
        const Platform platform;
        Feed(FeedSpecification const * specification) :
            specification(specification),
            direction(specification->direction),
            url(specification->url),
            platform(specification->platform),
            capacity(specification->capacity),
            resolution(specification->resolution),
            phred_offset(specification->phred_offset),
            end_of_file(false) {
            hfile = specification->hfile;
        };
        virtual ~Feed() {
        };
        virtual void join() = 0;
        virtual void start() = 0;
        virtual void stop() = 0;
        virtual void open() = 0;
        virtual void close() = 0;
        virtual bool pull(Segment& segment) = 0;
        virtual void push(Segment& segment) = 0;
        virtual bool peek(Segment& segment, const size_t& position) = 0;
        virtual inline void flush() = 0;
        virtual inline void replenish() = 0;
        virtual void calibrate(FeedSpecification const * specification) = 0;
        virtual unique_lock<mutex> acquire_pull_lock() = 0;
        virtual unique_lock<mutex> acquire_push_lock() = 0;
        virtual inline bool opened() = 0;
        virtual inline bool exhausted() = 0;
        void set_thread_pool(htsThreadPool* pool) {
            thread_pool = pool;
        };
        const size_t& index() const {
            return specification->index;
        };

    protected:
        hFILE* hfile;
        htsThreadPool* thread_pool;
        size_t capacity;
        size_t resolution;
        uint8_t phred_offset;
        bool end_of_file;
};
template <class T> class CyclicBuffer {
    template <typename U> friend ostream& operator<<(ostream& o, const CyclicBuffer<U>& buffer);

    public:
        CyclicBuffer (
            const IoDirection& direction,
            const size_t& capacity,
            const size_t& resolution) :

            _direction(direction),
            _capacity(0),
            _resolution(0),
            _next(-1),
            _vacant(0) {

            calibrate(capacity, resolution);
        };
        virtual ~CyclicBuffer() {
        };
        void increment() {
            // add_cache_record
            // if the cache is empty
            // next record is the vacant record we just populated
            if(_next < 0)
                _next = _vacant;

            // increment vacancy pointer on a circular buffer
            _vacant = (_vacant + 1) % _capacity;

            // if the vacancy pointer is the next record there is no more space
            // vacant becomes -1, or unavailable
            if(_vacant == _next)
                _vacant = -1;
        };
        void decrement() {
            // evacuate_cache_record
            // if there is no vacant record, next becomes the first vacant record
            if(_vacant < 0)
                _vacant = _next;

            // next record iterates on a circular buffer
            _next = (_next + 1) % _capacity;

            // if the next is the vacant record that means we have nothing in the buffer
            // and so next record is -1, i.e. undefined.
            if(_next == _vacant)
                _next = -1;
        };
        void increment_front() {
            // encode the segment in the front of the buffer
            if(_next < 0) {
                _next = _vacant;
                _vacant = (_vacant + 1) % _capacity;
            } else {
                _next = (_next - 1) % _capacity;
            }

            if(_vacant == _next)
                _vacant = -1;
        };
        inline T* vacant() const {
            return cache[_vacant];
        };
        inline T* next() const {
            return cache[_next];
        };
        inline size_t size() const {
            if(_next < 0) return 0;
            if(_vacant < 0) return _capacity;
            return (_vacant - _next) % _capacity;
        };
        inline size_t available() const {
            if(_next < 0) return _capacity;
            if(_vacant < 0) return 0;
            return (_next - _vacant) % _capacity;
        };
        const inline size_t& capacity() const {
            return _capacity;
        };
        const inline size_t& resolution() const {
            return _resolution;
        };
        const inline IoDirection& direction() const {
            return _direction;
        };
        inline bool is_full() const {
            return _vacant < 0;
        };
        inline bool is_not_full() const {
            return _vacant >= 0;
        };
        inline bool is_empty() const {
            return _next < 0;
        };
        inline bool is_not_empty() const {
            return _next >= 0;
        };
        void sync(CyclicBuffer<T>* other) {
            while(size() % resolution() > 0) {
                T* switching = other->cache[other->_next];
                other->cache[other->_next] = cache[_vacant];
                other->decrement();
                cache[_next] = switching;
                increment();
            }
        };
        virtual void calibrate(const size_t& capacity, const size_t& resolution);

    private:
        const IoDirection _direction;
        size_t _capacity;
        size_t _resolution;
        int _next;
        int _vacant;
        vector< T* > cache;
        int index;
};
template <class T> class BufferedFeed : public Feed {
    public:
        BufferedFeed(FeedSpecification const * specification) :
            Feed(specification),
            kbuffer({ 0, 0, NULL }),
            buffer(new CyclicBuffer<T>(direction, capacity, resolution)),
            queue(new CyclicBuffer<T>(direction, capacity, resolution)),
            started(false) {
            ks_terminate(kbuffer);
        };
        virtual ~BufferedFeed() {
            ks_free(kbuffer);
            delete queue;
            delete buffer;
        };
        void join() {
            feed_thread.join();
        };
        void start() {
            if(!started) {
                started = true;
                feed_thread = thread(&BufferedFeed::run, this);
            }
        };
        void stop() {
            lock_guard<mutex> feed_lock(queue_mutex);
            end_of_file = true;
            flushable.notify_one();
        };
        bool pull(Segment& segment) {
            if(queue->is_not_empty()) {
                decode(queue->next(), segment);
                queue->decrement();

                // queue is one element smaller, if its empty notify the replenish thread
                if(queue->is_empty()) {
                    replenishable.notify_one();
                }
                return true;
            }
            return false;
        };
        void push(Segment& segment) {
            encode(queue->vacant(), segment);
            queue->increment();

            // queue is one element bigger, if its full notify the flush thread
            if(is_ready_to_flush()) {
                flushable.notify_one();
            }
        };
        bool peek(Segment& segment, const size_t& position) {
            if(queue->size() > position) {
                decode(queue->next() + position, segment);
                return true;
            } else {
                segment.clear();
            }
            return false;
        };
        inline void flush() {
            unique_lock<mutex> buffer_lock(buffer_mutex);
            if(buffer->is_not_empty()) {
                // if(!opened()) open();
                empty_buffer();
            }

            // buffer is empty, wait for the queue to be full
            // and switch between buffer and queue
            unique_lock<mutex> queue_lock(queue_mutex);
            flushable.wait(queue_lock, [this](){ return is_ready_to_flush(); });

            // buffer is empty and queue is full
            // switch between buffer and queue
            switch_buffers();

            // now queue is empty and buffer is full
            // notify threads waiting for the queue to have available space
            queue_not_full.notify_all();
        };
        inline void replenish() {
            unique_lock<mutex> buffer_lock(buffer_mutex);
            if(buffer->is_not_full()) {
                // if(!opened()) open();
                fill_buffer();
            }

            if(buffer->is_not_empty()) {
                // the buffer is not empty, wait for the queue to be empty
                // and switch between the buffer and the queue
                unique_lock<mutex> queue_lock(queue_mutex);
                replenishable.wait(queue_lock, [this](){ return is_ready_to_replenish(); });

                // buffer is not empty and queue is empty
                // switch between buffer and queue
                switch_buffers();

                // now queue is not empty and buffer is empty
                // notify threads waiting for the queue to be available
                queue_not_empty.notify_all();
            }
        };
        void calibrate(FeedSpecification const * specification) {
            unique_lock<mutex> buffer_lock(buffer_mutex);
            unique_lock<mutex> queue_lock(queue_mutex);

            if(capacity != specification->capacity || resolution != specification->resolution) {
                if(specification->capacity > capacity) {
                    capacity = specification->capacity;
                    resolution = specification->resolution;
                    queue->calibrate(capacity, resolution);
                    buffer->calibrate(capacity, resolution);

                    /*  sync queue
                        move elements from buffer to queue until queue is aligned
                        with the new resolution
                    */
                    queue->sync(buffer);

                    /*  sync buffer
                        now make sure the buffer is filled which will align it
                    */
                    fill_buffer();
                } else {
                    throw InternalError("can not reduce buffer size");
                }
            }
        };
        unique_lock<mutex> acquire_pull_lock() {
            unique_lock<mutex> queue_lock(queue_mutex);
            queue_not_empty.wait(queue_lock, [this]() { return is_ready_to_pull(); });
            return queue_lock;
        };
        unique_lock<mutex> acquire_push_lock() {
            unique_lock<mutex> queue_lock(queue_mutex);
            queue_not_full.wait(queue_lock, [this]() { return is_ready_to_push(); });
            return queue_lock;
        };
        inline bool exhausted() {
            return queue->is_empty() && buffer->is_empty() && end_of_file;
        };

    protected:
        kstring_t kbuffer;
        CyclicBuffer<T>* buffer;
        CyclicBuffer<T>* queue;
        virtual inline void encode(T* record, const Segment& segment) const = 0;
        virtual inline void decode(const T* record, Segment& segment) = 0;
        virtual inline void fill_buffer() = 0;
        virtual inline void empty_buffer() = 0;
        inline void switch_buffers() {
            CyclicBuffer<T>* tmp = buffer;
            buffer = queue;
            queue = tmp;
        };
        inline bool is_ready_to_pull() {
            // should have enough to pull all segments
            return queue->is_not_empty() || (end_of_file && queue->is_empty() && buffer->is_empty());
        };
        inline bool is_ready_to_push() {
            // should have enough to push all segments
            return queue->is_not_full();
        };
        inline bool is_ready_to_flush() {
            return queue->is_full() || end_of_file;
        };
        inline bool is_ready_to_replenish() {
            return queue->is_empty();
        };

    private:
        bool started;
        thread feed_thread;
        mutex buffer_mutex;
        mutex queue_mutex;
        condition_variable queue_not_empty;
        condition_variable replenishable;
        condition_variable queue_not_full;
        condition_variable flushable;
        void run() {
            switch(direction) {
                case IoDirection::IN: {
                    do {
                        replenish();
                    } while(!end_of_file);
                    break;
                };
                case IoDirection::OUT: {
                    while(!exhausted()) {
                        flush();
                    }
                    break;
                };
            }
            close();
        };
};
class FastqFeed : public BufferedFeed<FastqRecord> {
    friend class Channel;

    public:
        FastqFeed(FeedSpecification const * specification) :
            BufferedFeed<FastqRecord>(specification),
            bgzf_file(NULL) {
        };
        void open() {
            if(!opened()) {
                // Enable multi-threading (when compiled with -DBGZF_MT) via a shared
                // thread pool.  This means both encoder and decoder can balance
                // usage across a single pool of worker jobs.
                // @param fp          BGZF file handler; must be opened for writing
                // @param pool        The thread pool (see hts_create_threads)
                // int bgzf_thread_pool(BGZF *fp, struct hts_tpool *pool, int qsize);

                switch(direction) {
                    case IoDirection::IN: {
                        bgzf_file = bgzf_hopen(hfile, "r");
                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + url + " for reading");
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        if(!strcmp(specification->url.compression(), "gz")) {
                            bgzf_file = bgzf_hopen(hfile, "wg");
                        } else {
                            bgzf_file = bgzf_hopen(hfile, "wu");
                        }
                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + url + " for writing");
                        }
                        break;
                    };
                }
            }
        };
        void close() {
            if(opened()) {
                bgzf_close(bgzf_file);
                bgzf_file = NULL;

                kseq_destroy(kseq);
                kseq = NULL;
            }
        };
        inline bool opened() {
            return bgzf_file != NULL;
        };

    protected:
        BGZF* bgzf_file;
        kseq_t* kseq;
        inline void encode(FastqRecord* record, const Segment& segment) const {
            record->decode(segment);
        };
        inline void decode(const FastqRecord* record, Segment& segment) {
            record->encode(segment);
        };
        inline void fill_buffer() {
            while(buffer->is_not_full()) {
             /* >=0  length of the sequence (normal)
                -1   end-of-file
                -2   truncated quality string */
                if(kseq_read(kseq) < 0) {
                    end_of_file = true;
                    break;
                } else {
                    buffer->vacant()->decode(kseq, phred_offset);
                    buffer->increment();
                }
            }
        };
        inline void empty_buffer() {
            ks_clear(kbuffer);
            while(buffer->is_not_empty()) {
                // encode FastqRecord records to text
                FastqRecord* record = buffer->next();
                record->encode(kbuffer, phred_offset);
                buffer->decrement();
            }

            // write the text encoded fastq records to the stream
            if(bgzf_write(bgzf_file, kbuffer.s, kbuffer.l) < 0) {
                throw IOError("error writing to " + url);
            }
        };
};
class HtsFeed : public BufferedFeed<bam1_t> {
    friend class Factory;
    friend class Channel;

    public:
        HtsFeed(FeedSpecification const * specification) :
            BufferedFeed<bam1_t>(specification),
            hts_file(NULL) {

            header.hd.set_alignment_sort_order(HtsSortOrder::UNKNOWN);
            header.hd.set_alignment_grouping(HtsGrouping::QUERY);
            for(const auto& record : specification->program_by_id) {
                header.add_program(record.second);
            }
            for(const auto& record : specification->read_group_by_id) {
                header.add_read_group(record.second);
            }
        };
        void open() {
            if(!opened()) {
                switch(direction) {
                    case IoDirection::IN: {
                        hts_file = hts_hopen(hfile, url.c_str(), "r");
                        if(hts_file != NULL) {
                            hts_set_thread_pool(hts_file, thread_pool);
                            header.decode(hts_file);
                        } else {
                            throw IOError("failed to open " + url + " for reading");
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        switch(url.type()) {
                            case FormatType::SAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "w");
                                if(hts_file) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::BAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "wb");
                                if(hts_file) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::CRAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "wc");
                                if(hts_file) {
                                    hts_file->format.version.major = 3;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            default:
                                break;
                        }
                        if (hts_file != NULL) {
                            hts_set_thread_pool(hts_file, thread_pool);
                            header.hd.set_version(&(hts_file->format));
                            header.assemble();
                            header.encode(hts_file);
                        } else {
                            throw IOError("failed to open " + url + " for writing");
                        }
                        break;
                    };
                }
            }
        };
        void close() {
            if(opened()) {
                hts_close(hts_file);
                hts_file = NULL;
            }
        };
        inline bool opened() {
            return hts_file != NULL;
        };
        HtsHeader& get_header() {
            return header;
        };

    protected:
        HtsHeader header;
        htsFile* hts_file;
        inline void encode(bam1_t* record, const Segment& segment) const {
            record->core.l_qname = segment.name.l + 1;
            record->core.l_qseq = segment.sequence.length;
            record->core.flag = segment.flag;

            size_t l_data = record->core.l_qname + ((record->core.l_qseq + 1)>>1) + record->core.l_qseq;
            record->l_data = l_data;

            // increase allocated space if necessary
            if(record->m_data < l_data) {
                record->m_data = l_data;
                kroundup32(record->m_data);
                record->data = (uint8_t*)realloc(record->data, record->m_data);
            }

            // write identifier to data block
            memcpy(bam_get_qname(record), segment.name.s, record->core.l_qname);

            // encode nucleotide byte BAM numeric encoding into nybble BAM numeric encoding
            uint8_t *bam_seq = bam_get_seq(record);
            for (int i = 0; i < record->core.l_qseq; i++) {
                bam1_seq_seti(bam_seq, i, segment.sequence.code[i]);
            }

            // encode the quality sequence
            memcpy(bam_get_qual(record), segment.sequence.quality, record->core.l_qseq);

            // encode the auxiliary fields into the record
            segment.auxiliary.encode(record);
        };
        inline void decode(const bam1_t* record, Segment& segment) {
            // decode identifier and write to segment
            kputsn(bam_get_qname(record), record->core.l_qname, &segment.name);

            // decode nucleotide sequence to buffer
            ks_clear(kbuffer);
            ks_resize(&kbuffer, record->core.l_qseq + 1);

            uint8_t* position = bam_get_seq(record);
            for(int32_t i = 0; i < record->core.l_qseq; i++) {
                // convert from 4bit BAM numeric encoding to 8bit BAM numeric encoding
                kbuffer.s[i] = bam_seqi(position, i);
                // kputc(bam_seqi(position, i), &kbuffer);
            }
            kbuffer.l = record->core.l_qseq;
            kbuffer.s[kbuffer.l] = '\0';
            segment.sequence.fill((uint8_t*)kbuffer.s, bam_get_qual(record), record->core.l_qseq);

            // decode flag
            segment.flag = record->core.flag;

            // decode auxiliary fields
            segment.auxiliary.decode(record);
        };
        inline void fill_buffer() {
            while(buffer->is_not_full()) {
                if(sam_read1(hts_file, header.hdr, buffer->vacant()) < 0) {
                    end_of_file = true;
                    close();
                    break;
                } else {
                    buffer->increment();
                }
            }
        };
        inline void empty_buffer() {
            while(buffer->is_not_empty()) {
                if(sam_write1(hts_file, header.hdr, buffer->next()) < 0) {
                    throw IOError("error writing to " + url);
                }
                buffer->decrement();
            }
        };
};
#endif /* PHENIQS_FEED_H */
