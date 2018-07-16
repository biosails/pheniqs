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

#include "environment.h"

Job* Environment::pop_from_queue() {
    if(!job_queue.empty()) {
        Job* job(job_queue.front());
        job_queue.pop_front();
        return job;
    } else { return NULL;}
};
void Environment::push_to_queue(Document operation) {
    if(operation.IsObject()) {
        string implementation(decode_value_by_key< string >("implementation", operation));
        Job* job(NULL);
        if(implementation == "demultiplex") {
            job = new Demultiplex(operation);
        } else {
            job = new Job(operation);
        }
        job->assemble();
        job_queue.emplace_back(job);
    } else { throw ConfigurationError("Job operation element is not a dictionary"); }
};
void Environment::execute_job(Job* job) {
    if(job != NULL) {
        job->compile();
        if(job->is_validate_only()) {
            job->describe(cerr);

        } else if(job->is_lint_only()) {
            job->print_ontology(cout);

        } else {
            job->execute();
            job->print_report(cerr);

        }
    }
};
void Environment::execute() {
    if(is_help_only()) {
        print_help(cerr);

    } else if(is_version_only()) {
        print_version(cerr);

    } else {
        push_to_queue(interface.operation());

        Job* job(NULL);
        while((job = pop_from_queue()) != NULL) {
            execute_job(job);
            delete job;
        }
    }
};
