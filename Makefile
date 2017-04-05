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


MAJOR_REVISON := 0
MINOR_REVISON := 9

-include environment.mk

# Construct version variable
GIT_VERSION := $(shell git describe --abbrev=40 --always 2> /dev/null)
ifdef GIT_VERSION
	PACKAGE_VERSION := $(MAJOR_REVISON).$(MINOR_REVISON).$(GIT_VERSION)
else
	PACKAGE_VERSION := $(MAJOR_REVISON).$(MINOR_REVISON)
endif

LIBS = -lhts -lz
CC = clang++
CFLAGS = -c -std=c++11 -O3 -Wall -Wsign-compare 
LDFLAGS = $(LIBS)

# Debug flags
# CFLAGS = -O0 -std=c++11 -Wall -Wsign-compare -g -Wconversion

SOURCES = \
	constant.cpp \
	model.cpp \
	feed.cpp \
	environment.cpp \
	pipeline.cpp \
	pheniqs.cpp

OBJECTS = \
	constant.o \
	model.o \
	feed.o \
	environment.o \
	pipeline.o \
	pheniqs.o

EXECUTABLE = pheniqs

# OS detection
# PLATFORM := $(shell uname -s)
# ifeq ($(PLATFORM),Linux)
# 	# CFLAGS += 
# endif
# ifeq ($(PLATFORM),Darwin)
# 	# CFLAGS += -stdlib=libc++
# endif

ifdef HTSLIB
	CFLAGS += -I$(HTSLIB)/include
	LDFLAGS += -L$(HTSLIB)/lib
endif

ifdef RAPIDJSON
	CFLAGS += -I$(RAPIDJSON)/include
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

# Dependencies
# Regenerate modules when header files they import change
constant.o: \
	constant.h

model.o: \
	constant.h \
	error.h \
	model.h

feed.o: \
	constant.h \
	error.h \
	model.h \
	feed.h

environment.o: \
	constant.h \
	error.h \
	model.h \
	feed.h \
	interface.h \
	configuration.h \
	version.h \
	environment.h

pipeline.o: \
	feed.h \
	environment.h \
	pipeline.h

pheniqs.o: \
	pipeline.h
