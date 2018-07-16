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

#ifndef PHENIQS_DEMULTIPLEX_H
#define PHENIQS_DEMULTIPLEX_H

#include "include.h"
#include "pipeline.h"
#include "accumulate.h"
#include "fastq.h"
#include "hts.h"
#include "decoder.h"
#include "metric.h"

class Demultiplex;
class DemultiplexPivot;

class Demultiplex : public Job {
    Demultiplex(Demultiplex const &) = delete;
    void operator=(Demultiplex const &) = delete;

    public:
        Demultiplex(Document& operation) :
            Job(operation),
            decoder_projection_query("/projection/decoder"),
            barcode_projection_query("/projection/barcode"),
            decoder_repository_query("/decoder"),
            end_of_input(false),
            thread_pool({NULL, 0}){
        };
        ~Demultiplex() override;
        inline bool display_distance() const {
            return decode_value_by_key< bool >("display distance", ontology);
        };
        void load() override {
            validate_url_accessibility();
            load_thread_pool();
            load_input();
            load_output();
            load_pivot();
        };
        void start();
        void stop();
        void execute() override;
        void describe(ostream& o) const override;
        void populate_decoder(DiscreteDecoder< Channel >& decoder);
        bool pull(Read& read);
        void print_compiled(ostream& o) const override {
            Document compiled;
            compiled.CopyFrom(ontology, compiled.GetAllocator());
            compiled.RemoveMember("decoder");
            compiled.RemoveMember("configuration url");
            compiled.RemoveMember("full command");
            compiled.RemoveMember("input feed");
            compiled.RemoveMember("input feed by segment");
            compiled.RemoveMember("operation");
            compiled.RemoveMember("output feed");
            compiled.RemoveMember("program");
            print_json(compiled, o);
        };

    protected:
        const Pointer decoder_projection_query;
        const Pointer barcode_projection_query;
        const Pointer decoder_repository_query;

        void manipulate() override;
        void validate() override;

    private:
        bool end_of_input;
        htsThreadPool thread_pool;
        list< DemultiplexPivot > pivot_array;
        list< Feed* > input_feed_by_index;
        list< Feed* > output_feed_by_index;
        vector< Feed* > input_feed_by_segment;
        unordered_map< URL, Feed* > output_feed_by_url;
        void compile_PG();
        void compile_input();
        void detect_input();
        void compile_decoder(Value& value, int32_t& index, const Value& default_decoder, const Value& default_barcode);
        void compile_decoder_group(const Value::Ch* key);
        void compile_output();
        void compile_transformation(Value& value);
        void compile_codec(Value& value, const Value& default_decoder, const Value& default_barcode);
        void compile_decoder_transformation(Value& value);
        bool infer_PU(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined=false);
        bool infer_ID(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined=false);
        void pad_url_array_by_key(const Value::Ch* key, Value& container, const int32_t& cardinality);
        void cross_validate_io();
        void validate_decoder_group(const Value::Ch* key);
        void validate_decoder(Value& value);

        void validate_url_accessibility();
        void load_thread_pool();
        void load_input();
        void load_output();
        void load_pivot();
        void populate_channel(Channel& channel);
        void finalize();

        void print_global_instruction(ostream& o) const;
        void print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const;
        void print_codec_instruction(const Value& value, const bool& plural, ostream& o) const;
        void print_channel_instruction(const Value& value, ostream& o) const;
        void print_feed_instruction(const Value::Ch* key, ostream& o) const;
        void print_input_instruction(ostream& o) const;
        void print_template_instruction(ostream& o) const;
        void print_codec_template(const Value& value, ostream& o) const;
        void print_multiplex_instruction(ostream& o) const;
        void print_molecular_instruction(ostream& o) const;
        void print_cellular_instruction(ostream& o) const;
};

class DemultiplexPivot {
    DemultiplexPivot(DemultiplexPivot const &) = delete;
    void operator=(DemultiplexPivot const &) = delete;

    public:
        const int32_t index;
        const Platform platform;
        const int32_t leading_segment_index;
        const int32_t input_segment_cardinality;
        const int32_t output_segment_cardinality;
        Read input;
        Read output;
        DiscreteDecoder< Channel >* multiplex;
        vector< Decoder* > molecular;
        vector< Decoder* > cellular;
        InputAccumulator input_accumulator;
        OutputAccumulator output_accumulator;
        DemultiplexPivot(Demultiplex& job, const int32_t& index);
        void start() {
            pivot_thread = thread(&DemultiplexPivot::run, this);
        };
        void join() {
            pivot_thread.join();
        };
        inline void clear() {
            input.clear();
            output.clear();
        };

    private:
        Demultiplex& job;
        thread pivot_thread;
        const bool disable_quality_control;
        const TemplateRule template_rule;
        inline void transform() {
            template_rule.apply(input, output);
            multiplex->decode(input, output);
            for(auto& decoder : molecular) {
                decoder->decode(input, output);
            }
            for(auto& decoder : cellular) {
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
        void load_multiplex_decoder();
        void load_molecular_decoder();
        void load_cellular_decoder();
        void run();
};

#endif /* PHENIQS_DEMULTIPLEX_H */
