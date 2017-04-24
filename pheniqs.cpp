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

#include "constant.h"
#include "pipeline.h"

int main(int argc, char** argv) {
    try {
        Environment environment(argc, argv);
        switch(environment.state) {
            case ProgramState::OK: {
                if (environment.validate_only) {
                    environment.describe(cout);

                } else if (environment.lint_only) {
                    environment.print_configuration(cout);

                } else {
                    Pipeline pipeline(environment);
                    pipeline.execute();
                }
                break;
            };

            case ProgramState::HELP:
                environment.print_help(cerr);
                break;

            case ProgramState::VERSION:
                environment.print_version(cerr);
                break;

            default:
                return static_cast< int >(environment.state);
                break;
        }

    } catch (InternalError error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::INTERNAL_ERROR);

    } catch (ConfigurationError error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::CONFIGURATION_ERROR);

    } catch (CommandLineError error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::COMMAND_LINE_ERROR);

    } catch (IOError error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::IO_ERROR);

    } catch (SequenceError error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::SEQUENCE_ERROR);

    } catch (exception error) {
        cerr << endl << error.what() << endl;
        return static_cast< int >(ProgramState::UNKNOWN_ERROR);
    }
    return static_cast< int >(ProgramState::OK);
};