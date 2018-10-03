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

Job::Job(Document& operation) try :
    operation(move(operation)),
    ontology(kObjectType),
    report(kObjectType),
    interactive(this->operation["interactive"]),
    schema_repository(this->operation["schema"]),
    projection_repository(this->operation["projection"]) {

    } catch(Error& error) {
        error.push("Job");
        throw;
};
void Job::clean() {
    clean_json_value(ontology, ontology);
    sort_json_value(ontology, ontology);
};
void Job::apply_default() {
    /* if the operation defines a default instruction overlay it on top of the ontology */
    Value::ConstMemberIterator reference = operation.FindMember("default");
    if(reference != operation.MemberEnd()) {
        merge_json_value(reference->value, ontology, ontology);
    }
};
void Job::apply_interactive() {
    overlay(interactive);
};
void Job::assemble() {
    /* if a URL to an instruction file was provided in the interactive instruction, load it and overlay on top of the ontology */
    URL configuration_url;
    if(decode_value_by_key< URL >("configuration url", configuration_url, interactive)) {
        configuration_url.normalize(IoDirection::IN);
        overlay(read_instruction_document(configuration_url));
    }
};
void Job::compile() {
    remove_disabled();
    clean();
    validate();
};
void Job::finalize() {
    if(decode_value_by_key< bool >("include compiled job", ontology)) {
        /*  add a copy of the job document */
        report.AddMember(
            Value("job", report.GetAllocator()).Move(),
            Value(ontology, report.GetAllocator()).Move(),
            report.GetAllocator()
        );
    }
};

void Job::print_ontology(ostream& o) const {
    print_json(ontology, o);
};
void Job::print_compiled(ostream& o) const {
    Document compiled;
    compiled.CopyFrom(ontology, compiled.GetAllocator());
    compiled.RemoveMember("application version");
    compiled.RemoveMember("program");
    compiled.RemoveMember("working directory");
    print_json(compiled, o, float_precision());
};
void Job::print_report() const {
    URL report_url(decode_value_by_key< URL >("report url", ontology));
    report_url.normalize(IoDirection::OUT);

    if(!report_url.is_dev_null()) {
        if(report_url.is_stdout()) {
            print_json(report, cout, float_precision());

        } else if(report_url.is_stderr()) {
            print_json(report, cerr, float_precision());

        } else {
            print_json(report, report_url.c_str(), float_precision());

        }
    }
};
void Job::describe(ostream& o) const {

};
void Job::overlay(const Value& instruction) {
    if(!instruction.IsNull()) {
        if(instruction.IsObject()) {
            if(!instruction.ObjectEmpty()) {
                Document merged;
                merged.CopyFrom(instruction, merged.GetAllocator());
                merge_json_value(ontology, merged, merged);
                ontology.Swap(merged);
            }
        } else { throw ConfigurationError("Job document root must be a dictionary"); }
    }
};
void Job::remove_disabled() {
    remove_disabled_from_json_value(ontology, ontology);
};
Document Job::read_instruction_document(const URL& url) {
    set< URL > visited;
    Document document(load_document_with_import(url, visited));
    return document;
};
Document Job::load_document_with_import(const URL& url, set< URL >& visited) {
    Document document(kNullType);
    if(url.is_readable()) {
        ifstream file(url.path());
        const string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
        file.close();
        if(!document.Parse(content.c_str()).HasParseError()) {
            const SchemaDocument* schema_document(get_schema_document("instruction:lax"));
            if(schema_document != NULL) {
                SchemaValidator validator(*schema_document);
                if(document.Accept(validator)) {
                    visited.emplace(url);
                    list< string > import;
                    if(decode_value_by_key< list< string > >("import", import, document)) {
                        Document aggregated(kNullType);
                        for(auto& record : import) {
                            URL import_url(expand_shell(record));

                            /* import url is resolved relative to the dirname of the importing document */
                            import_url.relocate_sibling(url);

                            /* To avoid cyclical import a url is only visited once,
                               the first time it is encountered on a depth first recursion.
                               TODO: This should really use the inode number and not the url */
                            if(!visited.count(import_url)) {
                                Document imported(load_document_with_import(import_url, visited));
                                merge_json_value(aggregated, imported, imported);
                                aggregated.Swap(imported);
                            }
                        }
                        merge_json_value(aggregated, document, document);
                    }
                    document.RemoveMember("import");

                } else {
                    Document error(encode_validation_error(validator, *find_schema("instruction:lax"), document));
                    encode_key_value("url", url, error, error);
                    throw ValidationError(error);
                }
            } else { throw InternalError("no schema instruction:lax"); }

        } else {
            string message(GetParseError_En(document.GetParseError()));
            message += " at position ";
            message += to_string(document.GetErrorOffset());
            throw ConfigurationError(message);
        }
    } else { throw ConfigurationError("unable to read job file from " + string(url)); }
    return document;
};
const Value* Job::find_projection(const string& key) const {
    const Value* element(NULL);
    if(projection_repository.IsObject()) {
        Value::ConstMemberIterator reference = projection_repository.FindMember(key.c_str());
        if(reference != projection_repository.MemberEnd()) {
            if(reference->value.IsObject()) {
                element = &reference->value;
            }
        }
    }
    return element;
};
const Value* Job::find_schema(const string& key) const {
    const Value* element(NULL);
    if(schema_repository.IsObject()) {
        Value::ConstMemberIterator reference = schema_repository.FindMember(key.c_str());
        if(reference != schema_repository.MemberEnd()) {
            if(reference->value.IsObject()) {
                element = &reference->value;
            }
        }
    }
    return element;
};
const SchemaDocument* Job::get_schema_document(const string& key) {
    const SchemaDocument* schema_document(NULL);
    auto record = schema_document_by_name.find(key);
    if(record != schema_document_by_name.end()) {
        schema_document = &record->second;
    } else {
        Value::ConstMemberIterator reference = schema_repository.FindMember(key.c_str());
        if(reference != schema_repository.MemberEnd()) {
            schema_document_by_name.emplace(make_pair(key, SchemaDocument(reference->value)));
            record = schema_document_by_name.find(key);
            if(record != schema_document_by_name.end()) {
                schema_document = &record->second;
            }
        }
    }
    return schema_document;
};
