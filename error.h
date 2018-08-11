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

/* Those are possible return values when pheniqs terminates */
enum class ErrorCode : int8_t {
    OK                         = 0,
    UNKNOWN_ERROR              = 1,
    INTERNAL_ERROR             = 2,
    CONFIGURATION_ERROR        = 3,
    OUT_OF_MEMORY_ERROR        = 4,
    COMMAND_LINE_ERROR         = 5,
    IO_ERROR                   = 6,
    SEQUENCE_ERROR             = 7,
    OVERFLOW_ERROR             = 8,
    CORRUPT_AUXILIARY_ERROR    = 9,
    JSON_VALIDATION_ERROR      = 10,
};

class Error : public exception {
    protected:
        Error(const string& name, const ErrorCode& code) :
            name(name),
            code(code) {
        };
        Error(const string& name, const ErrorCode& code, const string& message) :
            name(name),
            code(code),
            message(message) {
        };

    public:
        const string name;
        const ErrorCode code;
        list< string > stack;
        string message;
        void push(const string& where) {
            stack.emplace_back(where);
        };
        const char* what() const noexcept override {
            return name.c_str();
        };
        virtual ostream& describe(ostream& o) const {
            o << name;
            o << " : ";
            o << message;
            o << endl;
            return o;
        };
};

class OutOfMemoryError : public Error {
    public:
        OutOfMemoryError() :
            Error("Out of memory error", ErrorCode::OUT_OF_MEMORY_ERROR) {
        };
};

class InternalError : public Error {
    public:
        InternalError(const string& message) :
            Error("Internal error", ErrorCode::INTERNAL_ERROR, message) {
        };
};

class ConfigurationError : public Error {
    public:
        ConfigurationError(const string& message) :
            Error("Configuration error", ErrorCode::CONFIGURATION_ERROR, message) {
        };
};

class CommandLineError : public Error {
    public:
        CommandLineError(const string& message) :
            Error("Command line error", ErrorCode::COMMAND_LINE_ERROR, message) {
        };
};

class IOError : public Error {
    public:
        IOError(const string& message) :
            Error("IO error", ErrorCode::IO_ERROR, message) {
        };
};

class SequenceError : public Error {
    public:
        SequenceError(const string& message) :
            Error("Sequence error", ErrorCode::SEQUENCE_ERROR, message) {
        };
};

class OverflowError : public Error {
    public:
        OverflowError(const string& message) :
            Error("Overflow error", ErrorCode::OVERFLOW_ERROR, message) {
        };
};

class CorruptAuxiliaryError : public Error {
    public:
        CorruptAuxiliaryError(const string& message) :
            Error("Corrupt auxiliary error", ErrorCode::CORRUPT_AUXILIARY_ERROR, message) {
        };
};

#endif /* PHENIQS_ERROR_H */
