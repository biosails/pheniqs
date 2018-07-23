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
    report(kObjectType),
    projection_query("/projection") {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("Job :: " + error.message);

    } catch(exception& error) {
        throw InternalError("Job :: " + string(error.what()));
};
void Job::clean() {
    clean_json_value(ontology, ontology);
    sort_json_value(ontology, ontology);
};
void Job::assemble() {
    /* if the operation defines a default instruction overlay it on top of the ontology */
    Value::ConstMemberIterator reference = operation.FindMember("default");
    if(reference != operation.MemberEnd()) {
        overlay(reference->value);
    }

    reference = operation.FindMember("interactive");
    if(reference != operation.MemberEnd()) {
        const Value& interactive(reference->value);

        /* if a URL to an instruction file was provided in the interactive instruction, load it and overlay on top of the ontology */
        URL configuration_url;
        if(decode_value_by_key< URL >("configuration url", configuration_url, interactive)) {
            configuration_url.normalize(IoDirection::IN);
            overlay(read_instruction_document(configuration_url));
        }

        /* overlay the interactive instruction provided on the command line */
        overlay(interactive);
    }
    clean();
};
void Job::compile() {
    remove_disabled();
    manipulate();
    clean();
    validate();
};
void Job::print_ontology(ostream& o) const {
    print_json(ontology, o);
};
void Job::print_compiled(ostream& o) const {
    print_json(ontology, o);
};
void Job::print_report(ostream& o) const {
    print_json(report, o);
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
    remove_disabled_from_json_value(ontology);
};
Document Job::read_instruction_document(const URL& url) const {
    set< URL > visited;
    Document document(load_document_with_import(url, visited));
    return document;
};
Document Job::load_document_with_import(const URL& url, set< URL >& visited) const {
    Document document(kNullType);
    if(url.is_readable()) {
        ifstream file(url.path());
        const string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
        file.close();

        if(!document.Parse(content.c_str()).HasParseError()) {
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
    const Value* projection_dictionary(projection_query.Get(operation));
    if(projection_dictionary != NULL && projection_dictionary->IsObject()) {
        Value::ConstMemberIterator reference = projection_dictionary->FindMember(key.c_str());
        if(reference != ontology.MemberEnd()) {
            if(reference->value.IsObject()) {
                element = &reference->value;
            }
        }
    }
    return element;
};
