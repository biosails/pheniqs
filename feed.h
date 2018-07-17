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

#ifndef PHENIQS_FEED_H
#define PHENIQS_FEED_H

#include "include.h"
#include "proxy.h"
#include "read.h"
#include <htslib/thread_pool.h>

inline int align_to_resolution(const int& capacity, const int& resolution) {
    int aligned(static_cast< int >(capacity / resolution) * resolution);
    if(aligned < capacity) {
        aligned += resolution;
    }
    return aligned;
};

/* IO feed */
class Feed {
    public:
        const int32_t index;
        const URL url;
        const IoDirection direction;
        const uint8_t phred_offset;
        const Platform platform;
        Feed(const FeedProxy& proxy) :
            index(proxy.index),
            url(proxy.url),
            direction(proxy.direction),
            phred_offset(proxy.phred_offset),
            platform(proxy.platform),
            _capacity(proxy.capacity),
            _resolution(proxy.resolution),
            exhausted(false),
            hfile(proxy.hfile),
            thread_pool(NULL) {
        };
        virtual ~Feed() {
        };
        const inline int& capacity() const {
            return _capacity;
        };
        const inline int& resolution() const {
            return _resolution;
        };
        virtual void join() = 0;
        virtual void start() = 0;
        virtual void stop() = 0;
        virtual void open() = 0;
        virtual void close() = 0;
        virtual bool pull(Segment& segment) = 0;
        virtual void push(const Segment& segment) = 0;
        virtual bool peek(Segment& segment, const int& position) = 0;
        virtual inline bool flush() = 0;
        virtual inline bool replenish() = 0;
        virtual inline bool is_dev_null() {
            return url.is_dev_null();
        };
        virtual void calibrate_resolution(const int& resolution) = 0;
        virtual unique_lock< mutex > acquire_pull_lock() = 0;
        virtual unique_lock< mutex > acquire_push_lock() = 0;
        virtual inline bool opened() = 0;
        virtual void set_thread_pool(htsThreadPool* pool) {
            thread_pool = pool;
        };

    protected:
        int _capacity;
        int _resolution;
        bool exhausted;
        hFILE* hfile;
        htsThreadPool* thread_pool;
};

class NullFeed : public Feed {
    public:
        mutex null_mutex;
        NullFeed(const FeedProxy& proxy) :
            Feed(proxy) {
        };
        void join() override {
        };
        void start() override {
        };
        void stop() override {
        };
        void open() override {
        };
        void close() override {
        };
        bool pull(Segment& segment) override {
            return true;
        };
        void push(const Segment& segment) override {

        };
        bool peek(Segment& segment, const int& position) override {
            return false;
        };
        inline bool flush() override {
            return false;
        };
        inline bool replenish() override {
            return false;
        };
        void calibrate_resolution(const int& resolution) override {

        };
        unique_lock< mutex > acquire_pull_lock() override {
            unique_lock< mutex > queue_lock(null_mutex);
            return queue_lock;
        };
        unique_lock< mutex > acquire_push_lock() override {
            unique_lock< mutex > queue_lock(null_mutex);
            return queue_lock;
        };
        inline bool opened() override {
            return true;
        };
        void set_thread_pool(htsThreadPool* pool) override {

        };
};

template < class T > class CyclicBuffer {
    template < typename U > friend ostream& operator<<(ostream& o, const CyclicBuffer< U >& buffer);

    public:
        CyclicBuffer (
            const IoDirection& direction,
            const int& capacity,
            const int& resolution) :

            _direction(direction),
            _capacity(0),
            _resolution(resolution),
            _next(-1),
            _vacant(0) {

            calibrate_capacity(capacity);
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
        inline T* at(const int& position) const {
            if(position < size()) {
                return cache[(_next + position) % _capacity];
            } else {
                return NULL;
            }
        };
        inline int size() const {
            if(_next < 0) return 0;
            if(_vacant < 0) return _capacity;
            return (_vacant - _next) % _capacity;
        };
        inline int available() const {
            if(_next < 0) return _capacity;
            if(_vacant < 0) return 0;
            return (_next - _vacant) % _capacity;
        };
        const inline int& capacity() const {
            return _capacity;
        };
        const inline int& resolution() const {
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
        void sync(CyclicBuffer< T >* other) {
            while(size() % resolution() > 0) {
                T* switching = other->cache[other->_next];
                other->cache[other->_next] = cache[_vacant];
                other->decrement();
                cache[_next] = switching;
                increment();
            }
        };
        virtual int calibrate_capacity(const int& capacity);
        int calibrate_resolution(const int& resolution) {
            if(resolution != _resolution) {
                int aligned_capacity(align_to_resolution(_capacity, resolution));
                if(aligned_capacity > _capacity) {
                    calibrate_capacity(aligned_capacity);
                }
                _resolution = resolution;
            }
            return _capacity;
        };

    private:
        const IoDirection _direction;
        int _capacity;
        int _resolution;
        int _next;
        int _vacant;
        vector< T* > cache;
        int index;
};
template< typename T > ostream& operator<<(ostream& o, const CyclicBuffer< T >& buffer);

template < class T > class BufferedFeed : public Feed {
    private:
        inline void switch_buffer_and_queue() {
            CyclicBuffer< T >* tmp = buffer;
            buffer = queue;
            queue = tmp;
        };
        inline bool is_ready_to_flush() {
            return queue->is_full() || exhausted;
        };

    public:
        BufferedFeed(const FeedProxy& proxy) :
            Feed(proxy),
            kbuffer({ 0, 0, NULL }),
            buffer(new CyclicBuffer< T >(direction, proxy.capacity, proxy.resolution)),
            queue(new CyclicBuffer< T >(direction, proxy.capacity, proxy.resolution)),
            started(false) {
            ks_terminate(kbuffer);
        };
        virtual ~BufferedFeed() {
            ks_free(kbuffer);
            delete queue;
            delete buffer;
        };
        void join() override {
            feed_thread.join();
        };
        void start() override {
            if(!started) {
                started = true;
                feed_thread = thread(&BufferedFeed::run, this);
            }
        };
        void stop() override {
            lock_guard< mutex > feed_lock(queue_mutex);
            exhausted = true;
            flushable.notify_one();
        };
        bool pull(Segment& segment) override {
            /*  called in a safe context after acquire_pull_lock */
            if(queue->is_not_empty()) {
                decode(queue->next(), segment);
                queue->decrement();

                if(queue->is_empty()) {
                    /* wake up the replenishing thread */
                    replenishable.notify_one();
                }
                return true;
            }
            return false;
        };
        void push(const Segment& segment) override {
            encode(queue->vacant(), segment);
            queue->increment();

            if(is_ready_to_flush()) {
                flushable.notify_one();
            }
        };
        bool peek(Segment& segment, const int& position) override {
            if(queue->size() > position) {
                decode(queue->at(position), segment);
                return true;
            } else {
                segment.clear();
            }
            return false;
        };
        inline bool flush() override {
            /*  used by the consumer to empty the buffer into output */
            unique_lock< mutex > buffer_lock(buffer_mutex);
            flush_buffer();

            unique_lock< mutex > queue_lock(queue_mutex);
            flushable.wait(queue_lock, [this](){ return is_ready_to_flush(); });

            if(queue->is_not_empty()) {
                switch_buffer_and_queue();
                queue_not_full.notify_all();
                return true;
            } else {
                close();
                return false;
            }
        };
        inline bool replenish() override {
            /*  used by the producer to fill the buffer from the input */
            unique_lock< mutex > buffer_lock(buffer_mutex);
            replenish_buffer();

            unique_lock< mutex > queue_lock(queue_mutex);
            replenishable.wait(queue_lock, [this](){ return queue->is_empty(); });

            if(buffer->is_not_empty()) {
                switch_buffer_and_queue();
            } else {
                exhausted = true;
            }

            queue_not_empty.notify_all();
            return !exhausted;
        };
        void calibrate_resolution(const int& resolution) override {
            unique_lock< mutex > buffer_lock(buffer_mutex);
            unique_lock< mutex > queue_lock(queue_mutex);

            if(resolution != _resolution) {
                int aligned_capacity(align_to_resolution(_capacity, resolution));
                if(_capacity != aligned_capacity) {
                    queue->calibrate_resolution(resolution);
                    buffer->calibrate_resolution(resolution);
                    _capacity = aligned_capacity;
                    _resolution = resolution;

                    /*  sync queue
                        move elements from buffer to queue until
                        queue is aligned with the new resolution */
                    queue->sync(buffer);

                    /*  sync buffer
                        now make sure the buffer is filled which will align it */
                    replenish_buffer();

                } else { _resolution = resolution; }
            }
        };
        unique_lock< mutex > acquire_pull_lock() override {
            unique_lock< mutex > queue_lock(queue_mutex);
            queue_not_empty.wait(queue_lock, [this]() { return queue->is_not_empty() || exhausted; });
            return queue_lock;
        };
        unique_lock< mutex > acquire_push_lock() override {
            unique_lock< mutex > queue_lock(queue_mutex);
            queue_not_full.wait(queue_lock, [this]() { return queue->is_not_full(); });
            return queue_lock;
        };

    protected:
        kstring_t kbuffer;
        CyclicBuffer< T >* buffer;
        CyclicBuffer< T >* queue;
        virtual void encode(T* record, const Segment& segment) const = 0;
        virtual void decode(const T* record, Segment& segment) = 0;
        virtual void replenish_buffer() = 0;
        virtual void flush_buffer() = 0;

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
                    while(replenish());
                    break;
                };
                case IoDirection::OUT: {
                    while(flush());
                    break;
                };
                default:
                    break;
            }
        };
};
Value encode_value(const Feed& value, Document& document);
bool encode_key_value(const string& key, const list< Feed* >& value, Value& container, Document& document);
bool encode_key_value(const string& key, const vector< Feed* >& value, Value& container, Document& document);

#endif /* PHENIQS_FEED_H */
