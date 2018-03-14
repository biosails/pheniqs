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

#include "error.h"
#include "json.h"
#include "constant.h"
#include "url.h"
#include "sequence.h"
#include "auxiliary.h"
#include "model.h"
#include "segment.h"

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

static inline size_t align_capacity(const size_t& capacity, const size_t& resolution) {
    size_t aligned = size_t(capacity / resolution) * resolution;
    if(aligned < capacity) {
        aligned += resolution;
    }
    return aligned;
};

/* IO feed
*/
class Feed {
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
            /*  if the cache is empty then the next record 
                is the vacant record we just populated  */
            if(_next < 0) _next = _vacant;

            /*  increment vacancy pointer on the circular buffer  */
            _vacant = (_vacant + 1) % _capacity;

            /*  if the vacancy pointer is the next record there is 
                no more space and vacant becomes -1, or unavailable */
            if(_vacant == _next) _vacant = -1;
        };
        void decrement() {
            /*  if there is no vacant record, next becomes the first vacant record  */
            if(_vacant < 0) _vacant = _next;

            /*  iterate next on a circular buffer   */
            _next = (_next + 1) % _capacity;

            /*  if next is the vacant record that means the buffer is empty
                and so next becomes -1, or undefined    */
            if(_next == _vacant) _next = -1;
        };
        void increment_front() {
            /*  encode the segment in the front of the buffer
                This is not the usual policy and is used for calibration
            */
            if(_next < 0) {
                _next = _vacant;
                _vacant = (_vacant + 1) % _capacity;
            } else {
                _next = (_next - 1) % _capacity;
            }
            if(_vacant == _next) _vacant = -1;
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

                /*  queue is now one element smaller,
                    if its empty notify the replenishing thread */
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

            /*  queue is one element bigger,
                if its full notify the flushing thread */
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
                empty_buffer();
            }

            /*  buffer is empty, wait for the queue to be full
                and switch between buffer and queue     */
            unique_lock<mutex> queue_lock(queue_mutex);
            flushable.wait(queue_lock, [this](){ return is_ready_to_flush(); });

            /*  buffer is empty and queue is full
                switch between buffer and queue */
            switch_buffers();

            /*  now queue is empty and buffer is full
                notify threads waiting for the queue to have available space */
            queue_not_full.notify_all();
        };
        inline void replenish() {
            unique_lock<mutex> buffer_lock(buffer_mutex);
            if(buffer->is_not_full()) {
                fill_buffer();
            }

            if(buffer->is_not_empty()) {
                /*  the buffer is not empty, wait for the queue to be empty
                    and switch between the buffer and the queue */
                unique_lock<mutex> queue_lock(queue_mutex);
                replenishable.wait(queue_lock, [this](){ return is_ready_to_replenish(); });

                /*  buffer is not empty and queue is empty
                    switch between buffer and queue */
                switch_buffers();

                /*  now queue is not empty and buffer is empty
                    notify threads waiting for the queue to be available */
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
                        move elements from buffer to queue until 
                        queue is aligned with the new resolution */
                    queue->sync(buffer);

                    /*  sync buffer
                        now make sure the buffer is filled which will align it */
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
            return queue->is_not_empty() || (end_of_file && queue->is_empty() && buffer->is_empty());
        };
        inline bool is_ready_to_push() {
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
#endif /* PHENIQS_FEED_H */
