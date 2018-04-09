# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti <lior.galanti@nyu.edu>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


MAJOR_REVISON := 1
MINOR_REVISON := 1

# Construct version variable
GIT_VERSION := $(shell git describe --abbrev=40 --always 2> /dev/null)
ifdef GIT_VERSION
	PACKAGE_VERSION := $(MAJOR_REVISON).$(MINOR_REVISON).$(GIT_VERSION)
else
	PACKAGE_VERSION := $(MAJOR_REVISON).$(MINOR_REVISON)
endif

LIBS = -lhts -lz -lpthread
CC = clang++
CFLAGS = -c -std=c++11 -O3 -Wall -Wsign-compare 
LDFLAGS = $(LIBS)

SOURCES = \
	json.cpp \
	constant.cpp \
	url.cpp \
	interface.cpp \
	atom.cpp \
	transform.cpp \
	auxiliary.cpp \
	sequence.cpp \
	segment.cpp \
	specification.cpp \
	feed.cpp \
	fastq.cpp \
	hts.cpp \
	environment.cpp \
	accumulate.cpp \
	pipeline.cpp \
	pheniqs.cpp

OBJECTS = \
	json.o \
	constant.o \
	url.o \
	interface.o \
	atom.o \
	transform.o \
	auxiliary.o \
	sequence.o \
	segment.o \
	specification.o \
	feed.o \
	fastq.o \
	hts.o \
	environment.o \
	accumulate.o \
	pipeline.o \
	pheniqs.o

EXECUTABLE = pheniqs

ifdef PREFIX
	CFLAGS += -I$(PREFIX)/include
	LDFLAGS += -L$(PREFIX)/lib
endif

all: $(SOURCES) configuration.h version.h $(EXECUTABLE)

.cpp.o:
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# Regenerate version.h when PACKAGE_VERSION changes
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,clean-version))
	@echo Generate version.h with $(PACKAGE_VERSION)
	@echo '#define PHENIQS_VERSION "$(PACKAGE_VERSION)"' > $@

# Regenerate configuration.h when interface.json file changes
configuration.h: interface.json
	@echo Generate command line interface configuration
	@xxd -i interface.json > $@

clean-version:
	-@rm -f version.h

clean-configuration:
	-@rm -f configuration.h

clean: clean-version clean-configuration
	-@rm -f $(EXECUTABLE) $(OBJECTS)

install: pheniqs
	if( test ! -d $(PREFIX)/bin ) ; then mkdir -p $(PREFIX)/bin ; fi
	cp -f pheniqs $(PREFIX)/bin/pheniqs
	chmod a+x $(PREFIX)/bin/pheniqs

# Dependencies
# Regenerate modules when header files they import change
json.o: \
	error.h \
	json.h

constant.o: \
	error.h \
	json.h \
	constant.h

url.o: \
	error.h \
	json.h \
	constant.h \
	url.h

interface.o: \
	error.h \
	json.h \
	constant.h \
	interface.h

atom.o: \
	error.h \
	json.h \
	constant.h \
	atom.h

transform.o: \
	error.h \
	json.h \
	constant.h \
	transform.h

sequence.o: \
	error.h \
	json.h \
	constant.h \
	nucleotide.h \
	phred.h \
	transform.h \
	sequence.h

auxiliary.o: \
	error.h \
	json.h \
	constant.h \
	sequence.h \
	auxiliary.h

segment.o: \
	error.h \
	json.h \
	constant.h \
	nucleotide.h \
	sequence.h \
	auxiliary.h \
	segment.h

specification.o: \
	error.h \
	json.h \
	constant.h \
	url.h \
	atom.h \
	sequence.h \
	specification.h

accumulate.o: \
	error.h \
	json.h \
	constant.h \
	nucleotide.h \
	phred.h \
	sequence.h \
	segment.h \
	specification.h \
	accumulate.h

environment.o: \
	version.h \
	configuration.h \
	interface.h \
	error.h \
	json.h \
	constant.h \
	url.h \
	nucleotide.h \
	phred.h \
	atom.h \
	specification.h \
	environment.h

feed.o: \
	error.h \
	json.h \
	constant.h \
	url.h \
	sequence.h \
	specification.h \
	feed.h

fastq.o: \
	error.h \
	json.h \
	constant.h \
	nucleotide.h \
	phred.h \
	sequence.h \
	segment.h \
	specification.h \
	feed.h \
	fastq.h

hts.o: \
	error.h \
	json.h \
	constant.h \
	nucleotide.h \
	phred.h \
	atom.h \
	sequence.h \
	segment.h \
	specification.h \
	feed.h \
	hts.h

pipeline.o: \
	error.h \
	json.h \
	constant.h \
	url.h \
	nucleotide.h \
	phred.h \
	atom.h \
	accumulate.h \
	environment.h \
	feed.h \
	fastq.h \
	hts.h \
	pipeline.h

pheniqs.o: \
	error.h \
	json.h \
	constant.h \
	environment.h \
	pipeline.h