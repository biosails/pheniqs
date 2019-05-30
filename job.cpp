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

#include "job.h"

Job::Job(Document& operation) try :
    operation(move(operation)),
    interactive(this->operation["interactive"]),
    schema_repository(this->operation["schema"]),
    projection_repository(this->operation["projection"]),
    instruction(kObjectType),
    ontology(kObjectType),
    report(kObjectType) {

    } catch(Error& error) {
        error.push("Job");
        throw;
};
Job::~Job() {

};
void Job::assemble() {
    /* if a URL to an instruction file was provided in the interactive instruction,
       load it and overlay on top of the ontology */
    URL configuration_url;
    if(decode_value_by_key< URL >("configuration url", configuration_url, interactive)) {
        string buffer(configuration_url);
        expand_shell(buffer);
        normalize_standard_stream(buffer, IoDirection::IN);
        URL resolved(buffer);
        overlay_json_object(instruction, read_instruction_document(resolved));
        sort_json_value(instruction, instruction);
    }
};
void Job::compile() {
    ontology.CopyFrom(instruction, ontology.GetAllocator());
    remove_disabled_from_json_value(ontology, ontology);
    clean_json_object(ontology, ontology);
    validate();
};
void Job::execute() {
    load();
    start();
    stop();
    finalize();
};
void Job::write_report() const {
    URL report_url(decode_value_by_key< URL >("report url", ontology));
    if(!report_url.is_dev_null()) {
        if(report_url.is_stdout()) {
            print_json(report, cout, float_precision());

        } else if(report_url.is_stderr()) {
            print_json(report, cerr, float_precision());

        } else {
            print_json(report, report_url.path().c_str(), float_precision());

        }
    }
};
void Job::describe(ostream& o) const {

};
void Job::write_static_instruction(ostream& o) const {
    Document assembled;
    assembled.CopyFrom(instruction, assembled.GetAllocator());
    apply_interactive_ontology(assembled);
    sort_json_value(assembled, assembled);
    clean_json_object(assembled, assembled);
    print_json(assembled, o);
};
void Job::write_compiled_instruction(ostream& o) const {
    Document compiled;
    compiled.CopyFrom(ontology, compiled.GetAllocator());
    compiled.RemoveMember("application version");
    compiled.RemoveMember("program");
    compiled.RemoveMember("working directory");
    sort_json_value(compiled, compiled);
    print_json(compiled, o, float_precision());
};
void Job::validate() {

};
void Job::load() {

};
void Job::start() {

};
void Job::stop() {

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
void Job::apply_default_ontology(Document& document) const {
    /* if the operation defines a default instruction overlay it on top of the ontology */
    Value::ConstMemberIterator reference = operation.FindMember("default");
    if(reference != operation.MemberEnd()) {
        merge_json_value(reference->value, document, document);
    }
};
void Job::apply_interactive_ontology(Document& document) const {
    Document adjusted;
    adjusted.CopyFrom(interactive, adjusted.GetAllocator());
    adjusted.RemoveMember("configuration url");
    adjusted.RemoveMember("static only");
    adjusted.RemoveMember("validate only");
    adjusted.RemoveMember("compile only");
    overlay_json_object(document, adjusted);
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
Document Job::read_instruction_document(const URL& url) {
    set< URL > visited;
    Document document(load_instruction_with_import(url, visited));
    return document;
};
Document Job::load_instruction_with_import(const URL& url, set< URL >& visited) {
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
                                Document imported(load_instruction_with_import(import_url, visited));
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
    } else { throw ConfigurationError("unable to read instruction file from " + string(url.path())); }
    return document;
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
