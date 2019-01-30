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
        vector< T > element_by_index;

        RoutingClassifier(const Value& ontology) try :
            AccumulatingClassifier(decode_value_by_key< int32_t >("index", ontology)),
            decoded(NULL),
            unclassified(find_value_by_key("undetermined", ontology)),
            element_by_index(decode_value_by_key< vector< T > >("codec", ontology)) {

            decoded = &unclassified;

            } catch(Error& error) {
                error.push("RoutingClassifier");
                throw;
        };
        virtual inline void classify(const Read& input, Read& output) {
            ++(decoded->count);
            if(!output.qcfail()) {
                ++(decoded->pf_count);
            }
        };
        inline void finalize() override {
            for(auto& element : element_by_index) {
                this->classified_count += element.count;
                this->pf_classified_count += element.pf_count;
            }
            this->count = this->classified_count + unclassified.count;
            this->pf_count = this->pf_classified_count + unclassified.pf_count;

            for(auto& element : element_by_index) {
                element.finalize(*this);
            }
            unclassified.finalize(*this);
            AccumulatingClassifier::finalize();
        };
        RoutingClassifier< T >& operator+=(const RoutingClassifier< T >& rhs) {
            AccumulatingClassifier::operator+=(rhs);
            unclassified += rhs.unclassified;
            for(size_t index(0); index < element_by_index.size(); ++index) {
                element_by_index[index] += rhs.element_by_index[index];
            }
            return *this;
        };
        void encode(Value& container, Document& document) const override {
            AccumulatingClassifier::encode(container, document);

            Value unclassified_report(kObjectType);
            unclassified.encode(unclassified_report, document);
            container.AddMember("unclassified", unclassified_report.Move(), document.GetAllocator());

            if(!element_by_index.empty()) {
                Value element_report_array(kArrayType);
                for(auto& element : element_by_index) {
                    Value element_report(kObjectType);
                    element.encode(element_report, document);
                    element_report_array.PushBack(element_report.Move(), document.GetAllocator());
                }
                container.AddMember("classified", element_report_array.Move(), document.GetAllocator());
            }
        };
};

template < class T > class ReadGroupClassifier : public RoutingClassifier< T > {
    public:
        ReadGroupClassifier(const Value& ontology) try :
            RoutingClassifier< T >(ontology) {

            element_by_rg.reserve(this->element_by_index.size());
            for(auto& element : this->element_by_index) {
                element_by_rg.emplace(make_pair(string(element.rg.ID.s, element.rg.ID.l), &element));
            }

            } catch(Error& error) {
                error.push("ReadGroupClassifier");
                throw;
        };
        inline void classify(const Read& input, Read& output) override {
            this->decoded = &this->unclassified;
            if(ks_not_empty(input.RG())) {
                rg_id_buffer.assign(input.RG().s, input.RG().l);
                auto record = element_by_rg.find(rg_id_buffer);
                if(record != element_by_rg.end()) {
                    this->decoded = record->second;
                }
            }
        };

    protected:
        unordered_map< string, T* > element_by_rg;

    private:
        string rg_id_buffer;

};

#endif /* PHENIQS_CLASSIFY_H */
