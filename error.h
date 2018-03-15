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

#ifndef PHENIQS_ERROR_H
#define PHENIQS_ERROR_H

#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <exception>
#include <vector>
#include <string>

using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::setprecision;

class CommandLineError : public exception {
public:
    string message;
    CommandLineError(const string error) {
        message = error;
    };
    virtual const char* what() const throw() {
        return ("Command line error: " + message).c_str();
    };
};
class ConfigurationError : public exception {
public:
    string message;
    ConfigurationError(const string& error) {
        message.assign(error);
    };
    virtual const char* what() const throw() {
        return ("Configuration error: " + message).c_str();
    };
};
class SequenceError : public exception {
public:
    string message;
    SequenceError(const string& error) {
        message.assign(error);
    };
    virtual const char* what() const throw() {
        return ("Sequence error : " + message).c_str();
    };
};
class IOError : public exception {
public:
    string message;
    IOError(const string& error) {
        message.assign(error);
    };
    virtual const char* what() const throw() {
        return ("IO error : " + message).c_str();
    };
};
class InternalError : public exception {
public:
    string message;
    InternalError(const string& error) {
        message.assign(error);
    };
    virtual const char* what() const throw() {
        return ("Internal error : " + message).c_str();
    };
};

#endif /* PHENIQS_ERROR_H */