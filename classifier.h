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

template < class T > class Classifier : public AccumulatingClassifier {
    public:
        Classifier(const Value& ontology) try :
            AccumulatingClassifier(decode_value_by_key< int32_t >("index", ontology)) {

            } catch(Error& error) {
                error.push("Classifier");
                throw;
        };
        virtual ~Classifier() = default;
        virtual void classify(const Read& input, Read& output) = 0;
        virtual void finalize() = 0;
        // virtual Classifier< T >& operator+=(const Classifier< T >& rhs) = 0;
        virtual void encode(Value& container, Document& document) const = 0;
};

vector< string > decode_tag_ID_by_index(const Value& ontology);

template < class T > class RoutingClassifier : public AccumulatingClassifier {
    public:
        T* decoded;
        T unclassified;
        vector< T > tag_array;

        RoutingClassifier(const Value& ontology) try :
            AccumulatingClassifier(decode_value_by_key< int32_t >("index", ontology)),
            decoded(NULL),
            unclassified(ontology["undetermined"]),
            tag_array(decode_value_by_key< vector< T > >("codec", ontology)),
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
        inline void finalize() override {
            for(auto& element : tag_array) {
                this->classified_count += element.count;
                this->pf_classified_count += element.pf_count;
            }
            this->count = this->classified_count + unclassified.count;
            this->pf_count = this->pf_classified_count + unclassified.pf_count;

            for(auto& element : tag_array) {
                element.finalize(*this);
            }
            unclassified.finalize(*this);
            AccumulatingClassifier::finalize();
        };
        RoutingClassifier< T >& operator+=(const RoutingClassifier< T >& rhs) {
            AccumulatingClassifier::operator+=(rhs);
            unclassified += rhs.unclassified;
            for(size_t index(0); index < tag_array.size(); ++index) {
                tag_array[index] += rhs.tag_array[index];
            }
            return *this;
        };
        void encode(Value& container, Document& document) const override {
            AccumulatingClassifier::encode(container, document);

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
    protected:
        const bool multiplexing_classifier;
};

#endif /* PHENIQS_CLASSIFY_H */
