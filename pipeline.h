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

#ifndef PHENIQS_PIPELINE_H
#define PHENIQS_PIPELINE_H

#include "include.h"


#include "accumulate.h"
#include "environment.h"
#include "fastq.h"
#include "hts.h"
#include "decoder.h"

class Pivot;
class Pipeline;

class Pivot {
    Pivot(Pivot const &) = delete;
    void operator=(Pivot const &) = delete;

    public:
        const int32_t index;
        Read input;
        Read output;
        DiscreteDecoder< Channel >* multiplex;
        vector< Decoder* > molecular;
        vector< Decoder* > splitseq;
        InputAccumulator input_accumulator;
        OutputAccumulator output_accumulator;
        Pivot(Pipeline& pipeline, const Value& ontology, const int32_t& index);
        void start() {
            pivot_thread = thread(&Pivot::run, this);
        };
        void join() {
            pivot_thread.join();
        };
        inline void clear() {
            input.clear();
            output.clear();
        };

    private:
        Pipeline& pipeline;
        thread pivot_thread;
        const bool disable_quality_control;
        const TemplateRule template_rule;
        inline void transform() {
            template_rule.apply(input, output);
            multiplex->decode(input, output);
            for(auto& decoder : molecular) {
                decoder->decode(input, output);
            }
            for(auto& decoder : splitseq) {
                decoder->decode(input, output);
            }
        };
        inline void push() {
            multiplex->decoded->push(output);
        };
        inline void increment() {
            input_accumulator.increment(input);
            output_accumulator.increment(multiplex->decoded->index, output);
        };
        void load_multiplex_decoder(const Value& ontology);
        void load_molecular_decoder(const Value& ontology);
        void load_splitseq_decoder(const Value& ontology);
        void run();
};

class Pipeline {
    friend class Pivot;

    public:
        const Document& instruction;
        const Platform platform;
        const int32_t leading_segment_index;
        const bool disable_quality_control;
        const bool include_filtered;
        const int threads;
        const int32_t input_segment_cardinality;
        const int32_t output_segment_cardinality;
        Pipeline(const Document& instruction);
        virtual ~Pipeline();
        bool pull(Read& read);
        void execute();
        void populate_multiplex_decoder(DiscreteDecoder< Channel >& decoder);

    protected:
        void load();
        void start();
        void stop();
        void finalize();

    private:
        bool end_of_input;
        htsThreadPool thread_pool;
        list< Pivot > pivot_array;
        list< Feed* > input_feed_by_index;
        list< Feed* > output_feed_by_index;
        unordered_map< URL, Feed* > output_feed_by_url;
        vector< Feed* > input_feed_by_segment;
        InputAccumulator input_accumulator;
        OutputAccumulator output_accumulator;
        void validate_url_accessibility();
        void load_input();
        void load_output();
        void load_pivot();
        void populate_multiplex_channel(Channel& channel);
        void encode_report(ostream& o) const;
};

#endif /* PHENIQS_PIPELINE_H */
