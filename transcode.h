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

#ifndef PHENIQS_TRANSCODE_H
#define PHENIQS_TRANSCODE_H

#include "include.h"
#include "job.h"
#include "accumulator.h"
#include "fastq.h"
#include "hts.h"
#include "decoder.h"
#include "mdd.h"
#include "pamld.h"
#include "metric.h"
#include "multiplex.h"

class Transcode;
class TranscodePivot;

class CompoundClassifier {
    public:
        CompoundClassifier(const Value& ontology);
        vector< RoutingClassifier< Barcode >* > sample_classifier_array;
        vector< RoutingClassifier< Barcode >* > molecular_classifier_array;
        vector< RoutingClassifier< Barcode >* > cellular_classifier_array;
        void finalize();
        inline void classify(const Read& input, Read& output) {
            for(auto& classifier : sample_classifier_array) {
                classifier->classify(input, output);
            }
            for(auto& classifier : molecular_classifier_array) {
                classifier->classify(input, output);
            }
            for(auto& classifier : cellular_classifier_array) {
                classifier->classify(input, output);
            }
        };
        CompoundClassifier& operator+=(const CompoundClassifier& pivot);
        void encode(Value& container, Document& document) const;

    private:
        void load_multiplex_decoding(const Value& ontology);
        void load_sample_decoding(const Value& ontology);
        void load_sample_decoder(const Value& value);
        void load_molecular_decoding(const Value& ontology);
        void load_molecular_decoder(const Value& value);
        void load_cellular_decoding(const Value& ontology);
        void load_cellular_decoder(const Value& value);
};

class Transcode : public Job {
    friend class TranscodePivot;

    public:
        Transcode(Transcode const &) = delete;
        void operator=(Transcode const &) = delete;
        Transcode(Document& operation);
        ~Transcode() override;
        inline bool is_display_distance() const {
            return decode_value_by_key< bool >("display distance", ontology);
        };
        inline bool is_sense_input_layout() const {
            return decode_value_by_key< bool >("sense input layout", interactive);
        };
        bool pull(Read& read);
        void assemble() override;
        void compile() override;
        void describe() const override;

    protected:
        uint64_t count;
        uint64_t pf_count;
        double pf_fraction;
        void validate() override;
        void load() override;
        void start() override;
        void stop() override;
        void finalize() override;
        void apply_interactive_ontology(Document& document) const override;
        Transcode& operator+=(const TranscodePivot& pivot);

    private:
        bool end_of_input;
        int32_t decoded_nucleotide_cardinality;
        htsThreadPool thread_pool;
        list< TranscodePivot > pivot_array;
        list< Feed* > input_feed_by_index;
        list< Feed* > output_feed_by_index;
        vector< Feed* > input_feed_by_segment;
        unordered_map< URL, Feed* > output_feed_by_url;

        Multiplexer* multiplexer;
        vector< RoutingClassifier< Barcode >* > sample_classifier_array;
        vector< RoutingClassifier< Barcode >* > molecular_classifier_array;
        vector< RoutingClassifier< Barcode >* > cellular_classifier_array;

        void compile_PG();
        void compile_explicit_input();
        void compile_sensed_input();
        void compile_input();
        void apply_inheritance();
        void compile_barcode_decoding();

        void apply_topic_inheritance(const Value::Ch* key);
        void apply_decoder_inheritance(Value& value);
        void compile_topic(const Value::Ch* key);
        void compile_decoder(Value& value, int32_t& index, const Value& default_decoder, const Value& default_barcode);
        void compile_decoder_transformation(Value& value);
        void apply_repository_inheritence(const Value::Ch* key, Value& container, Document& document);
        void compile_output();
        void compile_template();
        void compile_transformation(Value& value);
        bool infer_PU(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined=false);
        bool infer_ID(const Value::Ch* key, string& buffer, Value& container, const bool& undetermined=false);
        void pad_url_array_by_key(const Value::Ch* key, Value& container, const int32_t& cardinality);
        void cross_validate_io();
        void compile_thread_model();
        void validate_decoder_group(const Value::Ch* key);
        void validate_decoder(Value& value);

        void validate_url_accessibility();
        void load_thread_pool();
        void load_decoding();
        void load_input();
        void load_output();
        void load_pivot();

        void load_multiplex_decoding();
        void load_sample_decoding();
        void load_sample_decoder(const Value& value);
        void load_molecular_decoding();
        void load_molecular_decoder(const Value& value);
        void load_cellular_decoding();
        void load_cellular_decoder(const Value& value);

        void print_global_instruction(ostream& o) const;
        void print_codec_group_instruction(const Value::Ch* key, const string& head, ostream& o) const;
        void print_codec_instruction(const Value& value, const bool& plural, ostream& o) const;
        void print_channel_instruction(const string& key, const Value& value, ostream& o) const;
        void print_input_instruction(ostream& o) const;
        void print_transform_instruction(ostream& o) const;
        void print_sample_instruction(ostream& o) const;
        void print_molecular_instruction(ostream& o) const;
        void print_cellular_instruction(ostream& o) const;
        void print_feed_instruction(const Value::Ch* key, ostream& o) const;
};

class TranscodePivot {
    public:
        TranscodePivot(TranscodePivot const &) = delete;
        void operator=(TranscodePivot const &) = delete;

        const int32_t index;
        const Platform platform;
        const int32_t leading_segment_index;
        const int32_t input_segment_cardinality;
        const int32_t output_segment_cardinality;
        Read input;
        Read output;
        Multiplexer multiplexer;
        vector< RoutingClassifier< Barcode >* > sample_classifier_array;
        vector< RoutingClassifier< Barcode >* > molecular_classifier_array;
        vector< RoutingClassifier< Barcode >* > cellular_classifier_array;
        TranscodePivot(Transcode& job, const int32_t& index);
        void start() {
            pivot_thread = thread(&TranscodePivot::run, this);
        };
        void join() {
            pivot_thread.join();
        };
        void finalize() {
            if(!sample_classifier_array.empty()) {
                for(auto& classifier : sample_classifier_array) {
                    classifier->finalize();
                }
            }
            if(!molecular_classifier_array.empty()) {
                for(auto& classifier : molecular_classifier_array) {
                    classifier->finalize();
                }
            }
            if(!cellular_classifier_array.empty()) {
                for(auto& classifier : cellular_classifier_array) {
                    classifier->finalize();
                }
            }
            multiplexer.finalize();
        };

    protected:
        void run() {
            while(job.pull(input)) {
                input.validate();
                if(!filter_incoming_qc_fail || !input.qcfail()) {
                    template_rule.apply(input, output);

                    /* decode sample, molecular and cellular barcodes */
                    for(auto& classifier : sample_classifier_array) {
                        classifier->classify(input, output);
                    }
                    for(auto& classifier : molecular_classifier_array) {
                        classifier->classify(input, output);
                    }
                    for(auto& classifier : cellular_classifier_array) {
                        classifier->classify(input, output);
                    }

                    /* complete output read assembly and push to output multiplexer */
                    output.flush();
                    multiplexer.push(output);
                }
                input.clear();
                output.clear();
            }
        };

    private:
        Transcode& job;
        thread pivot_thread;
        const bool filter_incoming_qc_fail;
        const TemplateRule template_rule;
        void load_decoding();
        void load_multiplex_decoding();
        void load_sample_decoding();
        void load_sample_decoder(const Value& value);
        void load_molecular_decoding();
        void load_molecular_decoder(const Value& value);
        void load_cellular_decoding();
        void load_cellular_decoder(const Value& value);
};

#endif /* PHENIQS_TRANSCODE_H */
