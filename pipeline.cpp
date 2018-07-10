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

#include "pipeline.h"

Job::Job(Document& node) :
    ontology(move(node)),
    report(kObjectType),
    operation_pointer("/operation"),
    operation_default_pointer("/operation/default"),
    operation_interactive_pointer("/operation/interactive"),
    operation_projection_pointer("/operation/projection") {

    Value* operation(operation_pointer.Get(ontology));
    if(operation != NULL) {
        if(operation->IsObject()) {
            string name(decode_value_by_key< string >("name", *operation));
        } else { throw ConfigurationError("Job operation element is not a dictionary"); }
    } else { throw ConfigurationError("Job ontology is missing an operation element"); }
};
void Job::compile() {
    compile_default();
    compile_from_url();
    compile_interactive();
    manipulate();
    clean();
    validate();
};
void Job::compile_default() {
    Value* operation_default(operation_default_pointer.Get(ontology));
    if(operation_default != NULL) {
        overlay(*operation_default);
    }
};
void Job::compile_from_url() {
    Value* operation_interactive(operation_interactive_pointer.Get(ontology));
    if(operation_interactive != NULL) {
        URL url;
        if(decode_value_by_key< URL >("configuration url", url, *operation_interactive)) {
            url.normalize(IoDirection::IN);
            overlay(load_document_from_url(url));
        }
    }
};
void Job::compile_interactive() {
    Value* operation_interactive(operation_interactive_pointer.Get(ontology));
    if(operation_interactive != NULL) {
        overlay(*operation_interactive);
    }
};
void Job::clean() {
    Value* operation(operation_pointer.Get(ontology));
    if(operation != NULL) {
        operation->RemoveMember("projection");
    }
    clean_json_value(ontology, ontology);
    sort_json_value(ontology, ontology);
};
void Job::overlay(const Value& value) {
    if(!value.IsNull()) {
        if(value.IsObject()) {
            if(!value.ObjectEmpty()) {
                Document merged;
                merged.CopyFrom(value, merged.GetAllocator());
                merge_json_value(ontology, merged, merged);
                ontology.Swap(merged);
                merged.SetNull();
            }
        } else { throw ConfigurationError("job element must be a dictionary"); }
    }
};
Document Job::load_document_from_url(const URL& url) {
    Document document(kNullType);
    if(url.is_readable()) {
        ifstream file(url.path());
        const string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
        file.close();
        if(!document.Parse(content.c_str()).HasParseError()) {
            apply_document_import(document);
        } else {
            string message(GetParseError_En(document.GetParseError()));
            message += " at position ";
            message += to_string(document.GetErrorOffset());
            throw ConfigurationError(message);
        }
    } else { throw ConfigurationError("unable to read job file from " + string(url)); }
    return document;
};
void Job::apply_document_import(Document& document) {
    list< string > import;
    if(decode_value_by_key< list< string > >("import", import, document)) {
        Document aggregated(kNullType);
        for(auto& record : import) {
            URL url(expand_shell(record));
            Document imported(load_document_from_url(url));
            merge_json_value(aggregated, imported, imported);
            aggregated.Swap(imported);
        }
        merge_json_value(aggregated, document, document);
    }
    document.RemoveMember("import");
};
