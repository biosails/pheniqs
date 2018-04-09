/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
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

#include <stdio.h>
#include <iostream>

#include "error.h"
#include "environment.h"
#include "pipeline.h"

using std::endl;
using std::cerr;
using std::cout;

int main(int argc, char** argv) {
    try {
        Environment environment(argc, (const char**)argv);

        if(environment.is_help_only()) {
            environment.print_help(cerr);

        } else if(environment.is_version_only()) {
            environment.print_version(cerr);

        } else if(environment.is_validate_only()) {
            environment.print_instruction_validation(cerr);

        } else if(environment.is_lint_only()) {
            environment.print_linted_instruction(cout);

        } else {
            Pipeline pipeline(environment);
            pipeline.execute();
        }

    } catch(InternalError& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::INTERNAL_ERROR);

    } catch(ConfigurationError& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::CONFIGURATION_ERROR);

    } catch(CommandLineError& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::COMMAND_LINE_ERROR);

    } catch(IOError& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::IO_ERROR);

    } catch(SequenceError& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::SEQUENCE_ERROR);

    } catch(exception& error) {
        cerr << error.what() << endl;
        return static_cast< int >(ProgramState::UNKNOWN_ERROR);

    }

    return static_cast< int >(ProgramState::OK);
};