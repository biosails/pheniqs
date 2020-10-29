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

#ifndef PHENIQS_CLASSIFY_H
#define PHENIQS_CLASSIFY_H

#include "include.h"
#include "selector.h"
#include "read.h"

enum class ClassifierType : int8_t {
    UNKNOWN     = -1,
    SAMPLE      =  0,
    CELLULAR    =  1,
    MOLECULAR   =  2,
};
string to_string(const ClassifierType& value);
bool from_string(const char* value, ClassifierType& result);
void to_kstring(const ClassifierType& value, kstring_t& result);
bool from_string(const string& value, ClassifierType& result);
ostream& operator<<(ostream& o, const ClassifierType& value);
template<> bool decode_value_by_key< ClassifierType >(const Value::Ch* key, ClassifierType& value, const Value& container);
template<> ClassifierType decode_value_by_key< ClassifierType >(const Value::Ch* key, const Value& container);

vector< string > decode_tag_id_by_index(const Value& ontology);

template < class T > class Classifier : public AccumulatingSelector {
    protected:
        T* decoded;
        T unclassified;
        vector< T > tag_array;
        const bool multiplexing_classifier;

    public:
        Classifier(const Value& ontology) try :
            AccumulatingSelector(decode_value_by_key< int32_t >("index", ontology)),
            decoded(NULL),
            unclassified(ontology["undetermined"]),
            tag_array(decode_value_by_key< vector< T > >("codec", ontology)),
            multiplexing_classifier(decode_value_by_key< bool >("multiplexing classifier", ontology)) {

            decoded = &unclassified;

            } catch(Error& error) {
                error.push("Classifier");
                throw;
        };
        Classifier(const Classifier< T >& other) :
            AccumulatingSelector(other),
            decoded(NULL),
            unclassified(other.unclassified),
            tag_array(other.tag_array),
            multiplexing_classifier(other.multiplexing_classifier) {

            decoded = &unclassified;
        };
        virtual inline void classify(const Read& input, Read& output) {
            ++(decoded->count);
            if(!output.qcfail()) {
                ++(decoded->pf_count);
            }
            if(multiplexing_classifier) {
                output.channel_index = decoded->index;
            }
        };
        virtual inline void collect(const Classifier& other) {
            AccumulatingSelector::collect(other);
            unclassified.collect(other.unclassified);
            for(size_t index(0); index < tag_array.size(); ++index) {
                tag_array[index].collect(other.tag_array[index]);
            }
        };
        inline void finalize() override {
            /* collect the counts from the tags to get the total */
            for(auto& element : tag_array) {
                this->classified_count += element.count;
                this->pf_classified_count += element.pf_count;
            }
            this->count = this->classified_count + unclassified.count;
            this->pf_count = this->pf_classified_count + unclassified.pf_count;

            /*  compute noise prior estimate
                first get the portion of reads that failed the noise filter from high confidence classified. */
            double estimated_noise_count(this->low_conditional_confidence_count);
            double confident_noise_ratio(estimated_noise_count / (estimated_noise_count + this->pf_classified_count));

            /* than assume that the same ratio applies to low confidnce reads */
            if(this->low_confidence_count > 0) {
                estimated_noise_count += double(this->low_confidence_count) * confident_noise_ratio;
            }
            this->estimated_noise_prior = estimated_noise_count / double(this->count);

            /* finalize the tags */
            double estimated_not_noise_prior(1.0 - this->estimated_noise_prior);
            for(auto& element : tag_array) {
                element.finalize(*this);
                element.estimated_concentration_prior = estimated_not_noise_prior * element.pf_pooled_classified_fraction;
            }
            unclassified.finalize(*this);

            /* finalize the selector */
            AccumulatingSelector::finalize();
        };
        void adjust_prior(Value& container, Document& document) {
            /* adjust the noise prior */
            encode_key_value("noise", estimated_noise_prior, container, document);

            /* build a map of barcode sequence to prior */
            unordered_map< string, double > concentration_prior_by_barcode(tag_array.size());
            kstring_t buffer({ 0, 0, NULL });
            for(auto& tag: tag_array) {
                tag.encode_iupac_ambiguity(buffer);
                concentration_prior_by_barcode.insert(make_pair<string, double>(string(buffer.s, buffer.l), double(tag.estimated_concentration_prior)));
                ks_clear(buffer);
            }
            ks_free(buffer);

            /* adjust the concentration prior on each barcode */
            Value::MemberIterator reference = container.FindMember("codec");
            if(reference != container.MemberEnd()) {
                list< string > barcode_segment;
                for(auto& record : reference->value.GetObject()) {
                    if(decode_value_by_key< list< string > >("barcode", barcode_segment, record.value)) {
                        string barcode_string;
                        for(auto& segment : barcode_segment) {
                            if(!barcode_string.empty()) {
                                barcode_string.append("-");
                            }
                            barcode_string.append(segment);

                            auto concentration_record = concentration_prior_by_barcode.find(barcode_string);
                            if(concentration_record != concentration_prior_by_barcode.end()) {
                                encode_key_value("concentration", concentration_record->second, record.value, document);
                            }
                        }
                    }
                }
            }
        };
        void encode(Value& container, Document& document) const override {
            AccumulatingSelector::encode(container, document);

            Value unclassified_report(kObjectType);
            unclassified.encode(unclassified_report, document);
            container.AddMember("unclassified", unclassified_report.Move(), document.GetAllocator());

            if(!tag_array.empty()) {
                Value element_report_array(kArrayType);
                for(auto& element : tag_array) {
                    Value element_report(kObjectType);
                    element.encode(element_report, document);
                    element_report_array.PushBack(element_report.Move(), document.GetAllocator());
                }
                container.AddMember("classified", element_report_array.Move(), document.GetAllocator());
            }
        };
};

#endif /* PHENIQS_CLASSIFY_H */
