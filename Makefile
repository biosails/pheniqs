# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2018  Lior Galanti
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

MAJOR_REVISON   := 2
MINOR_REVISON   := 0.6
PREFIX          := /usr/local
BIN_PREFIX      = $(PREFIX)/bin
INCLUDE_PREFIX  = $(PREFIX)/include
LIB_PREFIX      = $(PREFIX)/lib
ZSH_PREFIX      = $(PREFIX)/share/zsh

CC              = clang
CXX             = clang++
CPPFLAGS        += -Wall -Wsign-compare -Wdeprecated
CXXFLAGS        += -std=c++11 -O3
# LDFLAGS       +=
LIBS            += -lhts -lz -lbz2 -llzma
STATIC_LIBS     += $(LIB_PREFIX)/libhts.a $(LIB_PREFIX)/libz.a $(LIB_PREFIX)/libbz2.a $(LIB_PREFIX)/liblzma.a

# PHENIQS_VERSION, written into version.h and reported when executing `pheniqs --version`,
# is taken from the git describe if present. Otherwise it falls back to $(MAJOR_REVISON).$(MINOR_REVISON).
# Providing PHENIQS_VERSION to make on the command line during compilation will override both.
#
# If provided on the command line when executing make, the following dependecy variables will also be reported
# by `pheniqs --version`. This is mostly useful when building pheniqs statically so the user can tell which version
# of each dependecy is built into the finaly binary.
#
# PHENIQS_ZLIB_VERSION
# PHENIQS_BZIP2_VERSION
# PHENIQS_XZ_VERSION
# PHENIQS_LIBDEFLATE_VERSION
# PHENIQS_RAPIDJSON_VERSION
# PHENIQS_HTSLIB_VERSION
ifndef PHENIQS_VERSION
    PHENIQS_VERSION := $(shell [ -d .git ] && git describe --abbrev=40 --tags 2> /dev/null)
endif
ifndef PHENIQS_VERSION
    PHENIQS_VERSION := $(MAJOR_REVISON).$(MINOR_REVISON)
endif

PLATFORM := $(shell uname -s)

PHENIQS_SOURCES = \
	selector.cpp \
	atom.cpp \
	auxiliary.cpp \
	barcode.cpp \
	multiplex.cpp \
	decoder.cpp \
	classifier.cpp \
	naive.cpp \
	mdd.cpp \
	pamld.cpp \
	pipeline.cpp \
	fastq.cpp \
	feed.cpp \
	hts.cpp \
	interface.cpp \
	json.cpp \
	pheniqs.cpp \
	job.cpp \
	transcode.cpp \
	phred.cpp \
	proxy.cpp \
	read.cpp \
	sequence.cpp \
	transform.cpp \
	url.cpp

PHENIQS_OBJECTS = \
	selector.o \
	atom.o \
	auxiliary.o \
	barcode.o \
	multiplex.o \
	decoder.o \
	classifier.o \
	naive.o \
	mdd.o \
	pamld.o \
	pipeline.o \
	fastq.o \
	feed.o \
	hts.o \
	interface.o \
	json.o \
	pheniqs.o \
	job.o \
	transcode.o \
	phred.o \
	proxy.o \
	read.o \
	sequence.o \
	transform.o \
	url.o

PHENIQS_EXECUTABLE = pheniqs

ifdef PREFIX
    CPPFLAGS += -I$(INCLUDE_PREFIX)
    LDFLAGS += -L$(LIB_PREFIX)
endif

with-static = 0
with-libdeflate = 0
ifneq ('$(wildcard $(LIB_PREFIX)/libdeflate.a)','')
    ifeq ($(PLATFORM), Darwin)
        with-libdeflate = 1
    else ifeq ($(PLATFORM), Linux)
        ifneq ('$(wildcard $(LIB_PREFIX)/libdeflate.so)','')
            with-libdeflate = 1
        endif
    endif
endif

ifeq ($(with-libdeflate), 1)
    LIBS += -ldeflate
    STATIC_LIBS += $(LIB_PREFIX)/libdeflate.a
endif

ifeq ($(with-static), 1)
    LIBS = $(STATIC_LIBS)
endif


.PHONY: all
all: $(PHENIQS_SOURCES) generated $(PHENIQS_EXECUTABLE)

.PHONY: help
help:
	@printf '\nBuilding Pheniqs\n\
	\n\
	Code version : $(PHENIQS_VERSION)\n\
	\n\
	Make targets\n\
	\thelp      : Print this help.\n\
	\tall       : Build everything. This is the default if no target is specified.\n\
	\tclean     : Delete all generated and object files.\n\
	\tinstall   : Install pheniqs to $(PREFIX)\n\
	\tconfig    : Print the values of the influential variables and exit.\n\
	\ttest      : Run tests.\n\
	\t_pheniqs  : Generate the zsh completion script.\n\
	\n\
	Pheniqs depends on the following libraries:\n\
	\tzlib       : https://zlib.net\n\
	\tbz2        : http://www.bzip.org\n\
	\txz         : https://tukaani.org/xz\n\
	\tlibdeflate : https://github.com/ebiggers/libdeflate (optional)\n\
	\thtslib     : http://www.htslib.org\n\
	\trapidjson  : http://rapidjson.org\n\
	\n\
	libdeflate is a heavily optimized implementation of the DEFLATE algorithm.\n\
	Although pheniqs does not directly link to libdeflate, htslib can be optionaly linked to it when built,\n\
	which can significantly speed up reading and writing gzip compressed fastq files.\n\
	\n\
	To build pheniqs with a specific PREFIX, set the PREFIX variable when executing make.\n\
	Notice that you will need to specify it each time you execute make, not just when building.\n\
	For instance to build and install Pheniqs to /usr/local execute `make PREFIX=/usr/local && make install PREFIX=/usr/local`.\n\
	\n\
	To build pheniqs with a specific compiler, set the CXX variable when executing make.\n\
	For instance: `make CXX=/usr/local/bin/g++-7`.\n\
	\n\
	To build a statically linked binary set the `with-static` variable to 1.\n\
	This requires libhts.so, libz.so, libbz2.so, liblzma.so and optionally libdeflate.so \n\
	(libhts.a, libz.a, libbz2.a, liblzma.a and libdeflate.a on MacOS) to be available in LIB_PREFIX.\n\
	For instance: `make with-static=1`.\n\n'

.PHONY: config
config:
	$(if $(PHENIQS_VERSION),             $(info PHENIQS_VERSION             :  $(PHENIQS_VERSION)))
	$(if $(PHENIQS_VERSION),             $(info PHENIQS_VERSION             :  $(PHENIQS_VERSION)))
	$(if $(PLATFORM),                    $(info PLATFORM                    :  $(PLATFORM)))
	$(if $(PREFIX),                      $(info PREFIX                      :  $(PREFIX)))
	$(if $(BIN_PREFIX),                  $(info BIN_PREFIX                  :  $(BIN_PREFIX)))
	$(if $(INCLUDE_PREFIX),              $(info INCLUDE_PREFIX              :  $(INCLUDE_PREFIX)))
	$(if $(LIB_PREFIX),                  $(info LIB_PREFIX                  :  $(LIB_PREFIX)))
	$(if $(ZSH_PREFIX),                  $(info ZSH_PREFIX                  :  $(ZSH_PREFIX)))
	$(if $(CXX),                         $(info CXX                         :  $(CXX)))
	$(if $(CPPFLAGS),                    $(info CPPFLAGS                    :  $(CPPFLAGS)))
	$(if $(CXXFLAGS),                    $(info CXXFLAGS                    :  $(CXXFLAGS)))
	$(if $(LDFLAGS),                     $(info LDFLAGS                     :  $(LDFLAGS)))
	$(if $(LIBS),                        $(info LIBS                        :  $(LIBS)))
	$(if $(with-libdeflate),             $(info with-libdeflate             :  $(with-libdeflate)))
	$(if $(with-static),                 $(info with-static                 :  $(with-static)))
	$(if $(PHENIQS_ZLIB_VERSION),        $(info PHENIQS_ZLIB_VERSION        :  $(PHENIQS_ZLIB_VERSION)))
	$(if $(PHENIQS_BZIP2_VERSION),       $(info PHENIQS_BZIP2_VERSION       :  $(PHENIQS_BZIP2_VERSION)))
	$(if $(PHENIQS_XZ_VERSION),          $(info PHENIQS_XZ_VERSION          :  $(PHENIQS_XZ_VERSION)))
	$(if $(PHENIQS_LIBDEFLATE_VERSION),  $(info PHENIQS_LIBDEFLATE_VERSION  :  $(PHENIQS_LIBDEFLATE_VERSION)))
	$(if $(PHENIQS_RAPIDJSON_VERSION),   $(info PHENIQS_RAPIDJSON_VERSION   :  $(PHENIQS_RAPIDJSON_VERSION)))
	$(if $(PHENIQS_HTSLIB_VERSION),      $(info PHENIQS_HTSLIB_VERSION      :  $(PHENIQS_HTSLIB_VERSION)))

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

$(PHENIQS_EXECUTABLE): $(PHENIQS_OBJECTS)
	$(CXX) $(PHENIQS_OBJECTS) $(LDFLAGS) -pthread $(LIBS) -o $(PHENIQS_EXECUTABLE)

# Regenerate version.h when PHENIQS_VERSION changes
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PHENIQS_VERSION)",$(shell cat version.h)),,clean.version))
	$(info version.h generated with PHENIQS_VERSION $(PHENIQS_VERSION))
	$(if $(PHENIQS_VERSION),            @printf '#ifndef PHENIQS_VERSION_H\n'                                           >> $@)
	$(if $(PHENIQS_VERSION),            @printf '#define PHENIQS_VERSION_H\n\n'                                         >> $@)
	$(if $(PHENIQS_VERSION),            @printf '#define PHENIQS_VERSION "$(PHENIQS_VERSION)"\n'                        >> $@)
	$(if $(PHENIQS_ZLIB_VERSION),       @printf '#define PHENIQS_ZLIB_VERSION "$(PHENIQS_ZLIB_VERSION)"\n'              >> $@)
	$(if $(PHENIQS_BZIP2_VERSION),      @printf '#define PHENIQS_BZIP2_VERSION "$(PHENIQS_BZIP2_VERSION)"\n'            >> $@)
	$(if $(PHENIQS_XZ_VERSION),         @printf '#define PHENIQS_XZ_VERSION "$(PHENIQS_XZ_VERSION)"\n'                  >> $@)
	$(if $(PHENIQS_LIBDEFLATE_VERSION), @printf '#define PHENIQS_LIBDEFLATE_VERSION "$(PHENIQS_LIBDEFLATE_VERSION)"\n'  >> $@)
	$(if $(PHENIQS_RAPIDJSON_VERSION),  @printf '#define PHENIQS_RAPIDJSON_VERSION "$(PHENIQS_RAPIDJSON_VERSION)"\n'    >> $@)
	$(if $(PHENIQS_HTSLIB_VERSION),     @printf '#define PHENIQS_HTSLIB_VERSION "$(PHENIQS_HTSLIB_VERSION)"\n'          >> $@)
	$(if $(PHENIQS_VERSION),            @printf '\n#endif /* PHENIQS_VERSION_H */\n'                                    >> $@)

.PHONY: clean.version
clean.version:
	-@rm -f version.h

# Regenerate interface configuration.h from configuration.json
configuration.h: configuration.json
	$(and !$(shell ./tool/pheniqs-configuration-api.py header configuration.json > configuration.h), \
    $(info generating configuration.h) \
  )

# Regenerate zsh completion from configuration.json
_pheniqs: configuration.json
	$(and !$(shell ./tool/pheniqs-configuration-api.py zsh configuration.json > _pheniqs), \
    $(info zsh completion _pheniqs generated) \
  )

.PHONY: generated
generated: version.h configuration.h _pheniqs

.PHONY: clean.generated
clean.generated: clean.version
	-@rm -f configuration.h
	-@rm -f _pheniqs

.PHONY: clean.bin
clean.bin:
	-@rm -rf ./bin

.PHONY: clean.object
clean.object:
	-@rm -f $(PHENIQS_OBJECTS)

.PHONY: clean
clean: clean.generated clean.object clean.bin clean.test
	-@rm -f $(PHENIQS_EXECUTABLE)

.PHONY: install
install: pheniqs install.zsh_completion
	install pheniqs $(BIN_PREFIX)

.PHONY: install.zsh_completion
install.zsh_completion: _pheniqs
	$(and $(ZSH_PREFIX), \
    $(or $(wildcard $(ZSH_PREFIX)/site-functions),!$(shell mkdir -p $(ZSH_PREFIX)/site-functions)), \
    install _pheniqs $(ZSH_PREFIX)/site-functions \
)

.PHONY: uninstall.zsh_completion
uninstall.zsh_completion:
	-@rm -f $(ZSH_PREFIX)/site-functions/_pheniqs

.PHONY: uninstall
uninstall: uninstall.zsh_completion
	-@rm -f $(BIN_PREFIX)/pheniqs

.PHONY: test.api.configuration
test.api.configuration: all
	./test/api/configuration/run.sh

.PHONY: clean.test.api.configuration
clean.test.api.configuration:
	-@rm -rf test/api/configuration/result

.PHONY: test.pheniqs.BDGGG
test.pheniqs.BDGGG: all
	./test/BDGGG/run.sh

.PHONY: clean.test.pheniqs.BDGGG
clean.test.pheniqs.BDGGG:
	-@rm -rf test/BDGGG/result

.PHONY: test
test: test.api.configuration test.pheniqs.BDGGG

.PHONY: clean.test
clean.test: clean.test.api.configuration clean.test.pheniqs.BDGGG
	-@rm -rf test/api/configuration/result

# Dependencies

json.o: \
	error.h \
	kstring.h \
	json.h

url.o: \
	json.o \
	url.h

interface.o: \
	url.o \
	version.h \
	configuration.h \
	interface.h

atom.o: \
	json.o \
	atom.h

sequence.o: \
	json.o \
	phred.o \
	sequence.h

barcode.o: \
	sequence.o \
	selector.o \
	barcode.h

auxiliary.o: \
	atom.o \
	barcode.o \
	auxiliary.h

read.o: \
	auxiliary.o \
	read.h

selector.o: \
	json.o \
	selector.h

phred.o: \
	iupac.h \
	json.o \
	phred.h

proxy.o: \
	url.o \
	atom.o \
	proxy.h

feed.o: \
	proxy.o \
	read.o \
	feed.h

fastq.o: \
	feed.o \
	fastq.h

hts.o: \
	feed.o \
	hts.h

transform.o: \
	read.o \
	transform.h

multiplex.o: \
	feed.o \
	selector.o \
	multiplex.h

classifier.o: \
	selector.o \
	read.o \
	classifier.h

decoder.o: \
	classifier.o \
	transform.o \
	multiplex.o \
	decoder.h

naive.o: \
	decoder.o \
	naive.h

mdd.o: \
	decoder.o \
	mdd.h

pamld.o: \
	decoder.o \
	pamld.h

job.o: \
	json.o \
	url.o \
	job.h

transcode.o: \
	selector.o \
	fastq.o \
	hts.o \
	decoder.o \
	naive.o \
	mdd.o \
	pamld.o \
	metric.h \
	job.h \
	transcode.h

pipeline.o: \
	interface.o \
	transcode.o \
	pipeline.h

pheniqs.o: \
	pipeline.o
