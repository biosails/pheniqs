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

Job::Job(Document& operation) :
    operation(move(operation)),
    report(kObjectType),
    projection_query("/projection") {
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
void Job::remove_disabled() {
    remove_disabled_from_json_value(ontology);
};
void Job::compile() {
    remove_disabled();
    manipulate();
    clean();
    validate();
};
void Job::clean() {
    clean_json_value(ontology);
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
Document Job::read_instruction_document(const URL& url) const {
    Document document(kNullType);
    if(url.is_readable()) {
        ifstream file(url.path());
        const string content((istreambuf_iterator< char >(file)), istreambuf_iterator< char >());
        file.close();
        if(!document.Parse(content.c_str()).HasParseError()) {
            apply_instruction_import(document);
        } else {
            string message(GetParseError_En(document.GetParseError()));
            message += " at position ";
            message += to_string(document.GetErrorOffset());
            throw ConfigurationError(message);
        }
    } else { throw ConfigurationError("unable to read job file from " + string(url)); }
    return document;
};
void Job::apply_instruction_import(Document& instruction) const {
    list< string > import;
    if(decode_value_by_key< list< string > >("import", import, instruction)) {
        Document aggregated(kNullType);
        for(auto& record : import) {
            URL url(expand_shell(record));
            Document imported(read_instruction_document(url));
            merge_json_value(aggregated, imported, imported);
            aggregated.Swap(imported);
        }
        merge_json_value(aggregated, instruction, instruction);
    }
    instruction.RemoveMember("import");
};
