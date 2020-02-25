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
#include "accumulator.h"
#include "read.h"

template < class T > class RoutingClassifier : public AccumulatingClassifier {
    public:
        T* decoded;
        T unclassified;
        vector< T > tag_by_index;

        RoutingClassifier(const Value& ontology) try :
            AccumulatingClassifier(decode_value_by_key< int32_t >("index", ontology)),
            decoded(NULL),
            unclassified(find_value_by_key("undetermined", ontology)),
            tag_by_index(decode_value_by_key< vector< T > >("codec", ontology)),
            multiplexing_classifier(decode_value_by_key< bool >("multiplexing classifier", ontology)) {

            decoded = &unclassified;

            } catch(Error& error) {
                error.push("RoutingClassifier");
                throw;
        };
        RoutingClassifier(const RoutingClassifier< T >& other) :
            AccumulatingClassifier(other),
            decoded(NULL),
            unclassified(other.unclassified),
            tag_by_index(other.tag_by_index),
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
        inline void finalize() override {
            for(auto& element : tag_by_index) {
                this->classified_count += element.count;
                this->pf_classified_count += element.pf_count;
            }
            this->count = this->classified_count + unclassified.count;
            this->pf_count = this->pf_classified_count + unclassified.pf_count;

            for(auto& element : tag_by_index) {
                element.finalize(*this);
            }
            unclassified.finalize(*this);
            AccumulatingClassifier::finalize();
        };
        RoutingClassifier< T >& operator+=(const RoutingClassifier< T >& rhs) {
            AccumulatingClassifier::operator+=(rhs);
            unclassified += rhs.unclassified;
            for(size_t index(0); index < tag_by_index.size(); ++index) {
                tag_by_index[index] += rhs.tag_by_index[index];
            }
            return *this;
        };
        void encode(Value& container, Document& document) const override {
            AccumulatingClassifier::encode(container, document);

            Value unclassified_report(kObjectType);
            unclassified.encode(unclassified_report, document);
            container.AddMember("unclassified", unclassified_report.Move(), document.GetAllocator());

            if(!tag_by_index.empty()) {
                Value element_report_array(kArrayType);
                for(auto& element : tag_by_index) {
                    Value element_report(kObjectType);
                    element.encode(element_report, document);
                    element_report_array.PushBack(element_report.Move(), document.GetAllocator());
                }
                container.AddMember("classified", element_report_array.Move(), document.GetAllocator());
            }
        };
    protected:
        const bool multiplexing_classifier;
};

#endif /* PHENIQS_CLASSIFY_H */
