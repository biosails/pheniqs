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

#ifndef PHENIQS_ERROR_H
#define PHENIQS_ERROR_H

#include "include.h"

class InternalError : public exception {
    public:
        string message;
        InternalError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("Internal error : " + message).c_str();
        };
};

class OutOfMemoryError : public exception {
    public:
        const char* what() const throw() override {
            return "Out of memory error";
        };
};

class ConfigurationError : public exception {
    public:
        string message;
        ConfigurationError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("Configuration error: " + message).c_str();
        };
};

class CommandLineError : public exception {
    public:
        string message;
        CommandLineError(const string error) {
            message = error;
        };
        const char* what() const throw() override {
            return ("Command line error: " + message).c_str();
        };
};

class IOError : public exception {
    public:
        string message;
        IOError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("IO error : " + message).c_str();
        };
};

class SequenceError : public exception {
    public:
        string message;
        SequenceError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("Sequence error : " + message).c_str();
        };
};

class OverflowError : public exception {
    public:
        string message;
        OverflowError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("Identifiers error : " + message).c_str();
        };
};
class CorruptAuxiliaryError : public exception {
    public:
        string message;
        CorruptAuxiliaryError(const string& error) {
            message.assign(error);
        };
        const char* what() const throw() override {
            return ("Corrupt auxiliary error : " + message).c_str();
        };
};

#endif /* PHENIQS_ERROR_H */
