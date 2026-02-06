# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2018  Lior Galanti
# NYU Center for Genetics and System Biology
#
# Author: Lior Galanti <lior.galanti@nyu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.10)

# Project configuration
set(PHENIQS_MAJOR_REVISION 2)
set(PHENIQS_MINOR_REVISION 1.0)
set(PHENIQS_VERSION "${PHENIQS_MAJOR_REVISION}.${PHENIQS_MINOR_REVISION}")

project(pheniqs VERSION ${PHENIQS_VERSION} LANGUAGES CXX)

# Try to get version from git if available
find_package(Git QUIET)
if(Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --abbrev=40 --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE PHENIQS_GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(PHENIQS_GIT_VERSION)
        set(PHENIQS_VERSION ${PHENIQS_GIT_VERSION})
    endif()
endif()

# Installation directories
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")
set(BIN_PREFIX "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Binary installation directory")
set(INCLUDE_PREFIX "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Include installation directory")
set(LIB_PREFIX "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Library installation directory")
set(ZSH_PREFIX "${CMAKE_INSTALL_PREFIX}/share/zsh" CACHE PATH "Zsh completion installation directory")

# Compiler settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wsign-compare -Wdeprecated -O3")

# Use clang if available (matching Makefile default)
if(NOT CMAKE_CXX_COMPILER)
    find_program(CLANGXX_COMPILER clang++)
    if(CLANGXX_COMPILER)
        set(CMAKE_CXX_COMPILER ${CLANGXX_COMPILER})
    endif()
endif()

# Static build option
option(BUILD_STATIC "Build with static libraries" OFF)

# Find required packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Find htslib (no official CMake config in older versions, so we'll find it manually)
find_library(HTSLIB_LIBRARY NAMES hts PATHS ${LIB_PREFIX} NO_DEFAULT_PATH)
find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h PATHS ${INCLUDE_PREFIX} NO_DEFAULT_PATH)

if(HTSLIB_LIBRARY AND HTSLIB_INCLUDE_DIR)
    message(STATUS "Found htslib: ${HTSLIB_LIBRARY}")
else()
    # Try system paths as fallback
    find_library(HTSLIB_LIBRARY NAMES hts)
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h)
    if(HTSLIB_LIBRARY AND HTSLIB_INCLUDE_DIR)
        message(STATUS "Found htslib (system): ${HTSLIB_LIBRARY}")
    else()
        message(FATAL_ERROR "htslib not found. Please install htslib development files.")
    endif()
endif()

# Find zlib
find_package(ZLIB REQUIRED)

# Find bzip2
find_library(BZIP2_LIBRARY NAMES bz2 PATHS ${LIB_PREFIX} NO_DEFAULT_PATH)
find_path(BZIP2_INCLUDE_DIR NAMES bzlib.h PATHS ${INCLUDE_PREFIX} NO_DEFAULT_PATH)

if(NOT BZIP2_LIBRARY)
    find_package(BZip2 REQUIRED)
    set(BZIP2_LIBRARY BZip2::BZip2)
endif()

# Find liblzma
find_library(LZMA_LIBRARY NAMES lzma PATHS ${LIB_PREFIX} NO_DEFAULT_PATH)
find_path(LZMA_INCLUDE_DIR NAMES lzma.h PATHS ${INCLUDE_PREFIX} NO_DEFAULT_PATH)

if(NOT LZMA_LIBRARY)
    find_package(LibLZMA REQUIRED)
    set(LZMA_LIBRARY LibLZMA::LibLZMA)
endif()

# Find libdeflate (optional)
find_library(DEFLATE_LIBRARY NAMES deflate PATHS ${LIB_PREFIX} NO_DEFAULT_PATH)
find_path(DEFLATE_INCLUDE_DIR NAMES libdeflate.h PATHS ${INCLUDE_PREFIX} NO_DEFAULT_PATH)

if(DEFLATE_LIBRARY)
    message(STATUS "Found libdeflate: ${DEFLATE_LIBRARY}")
else()
    find_library(DEFLATE_LIBRARY NAMES deflate)
    if(DEFLATE_LIBRARY)
        message(STATUS "Found libdeflate (system): ${DEFLATE_LIBRARY}")
    else()
        message(STATUS "libdeflate not found (optional)")
        set(DEFLATE_LIBRARY "")
    endif()
endif()

# Collect libraries
set(PHENIQS_LIBRARIES 
    ${HTSLIB_LIBRARY}
    ZLIB::ZLIB
    ${BZIP2_LIBRARY}
    ${LZMA_LIBRARY}
    ${DEFLATE_LIBRARY}
)

# Handle static linking
if(BUILD_STATIC)
    # Find static versions
    find_library(HTSLIB_STATIC NAMES libhts.a PATHS ${LIB_PREFIX})
    find_library(ZLIB_STATIC NAMES libz.a PATHS ${LIB_PREFIX})
    find_library(BZIP2_STATIC NAMES libbz2.a PATHS ${LIB_PREFIX})
    find_library(LZMA_STATIC NAMES liblzma.a PATHS ${LIB_PREFIX})
    find_library(DEFLATE_STATIC NAMES libdeflate.a PATHS ${LIB_PREFIX})
    
    if(HTSLIB_STATIC AND ZLIB_STATIC AND BZIP2_STATIC AND LZMA_STATIC AND DEFLATE_STATIC)
        set(PHENIQS_LIBRARIES 
            ${HTSLIB_STATIC}
            ${ZLIB_STATIC}
            ${BZIP2_STATIC}
            ${LZMA_STATIC}
            ${DEFLATE_STATIC}
        )
        message(STATUS "Using static libraries")
    else()
        message(WARNING "Not all static libraries found, falling back to shared")
    endif()
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${INCLUDE_PREFIX})
include_directories(${HTSLIB_INCLUDE_DIR})

# Pheniqs source files
set(PHENIQS_SOURCES
    selector.cpp
    atom.cpp
    auxiliary.cpp
    barcode.cpp
    multiplex.cpp
    decoder.cpp
    classifier.cpp
    naive.cpp
    mdd.cpp
    pamld.cpp
    pipeline.cpp
    fastq.cpp
    feed.cpp
    hts.cpp
    interface.cpp
    json.cpp
    pheniqs.cpp
    job.cpp
    transcode.cpp
    phred.cpp
    proxy.cpp
    read.cpp
    sequence.cpp
    transform.cpp
    url.cpp
)

# Generate version.h
set(VERSION_HEADER "${CMAKE_BINARY_DIR}/version.h")
add_custom_command(
    OUTPUT "${VERSION_HEADER}"
    COMMAND ${CMAKE_COMMAND}
        -DPHENIQS_VERSION="${PHENIQS_VERSION}"
        -DOUTPUT_FILE="${VERSION_HEADER}"
        -DPHENIQS_ZLIB_VERSION="${PHENIQS_ZLIB_VERSION}"
        -DPHENIQS_BZIP2_VERSION="${PHENIQS_BZIP2_VERSION}"
        -DPHENIQS_XZ_VERSION="${PHENIQS_XZ_VERSION}"
        -DPHENIQS_LIBDEFLATE_VERSION="${PHENIQS_LIBDEFLATE_VERSION}"
        -DPHENIQS_RAPIDJSON_VERSION="${PHENIQS_RAPIDJSON_VERSION}"
        -DPHENIQS_HTSLIB_VERSION="${PHENIQS_HTSLIB_VERSION}"
        -P "${CMAKE_SOURCE_DIR}/cmake/GenerateVersionHeader.cmake"
    COMMENT "Generating version.h with version ${PHENIQS_VERSION}"
    BYPRODUCTS "${VERSION_HEADER}"
)

# Generate configuration.h from configuration.json
set(CONFIG_HEADER "${CMAKE_BINARY_DIR}/configuration.h")
add_custom_command(
    OUTPUT "${CONFIG_HEADER}"
    COMMAND python3 ${CMAKE_SOURCE_DIR}/tool/pheniqs-configuration-api.py header 
            ${CMAKE_SOURCE_DIR}/configuration.json > "${CONFIG_HEADER}"
    DEPENDS ${CMAKE_SOURCE_DIR}/configuration.json
    COMMENT "Generating configuration.h from configuration.json"
    BYPRODUCTS "${CONFIG_HEADER}"
)

# Generate zsh completion
set(ZSH_COMPLETION "${CMAKE_BINARY_DIR}/_pheniqs")
add_custom_command(
    OUTPUT "${ZSH_COMPLETION}"
    COMMAND python3 ${CMAKE_SOURCE_DIR}/tool/pheniqs-configuration-api.py zsh 
            ${CMAKE_SOURCE_DIR}/configuration.json > "${ZSH_COMPLETION}"
    DEPENDS ${CMAKE_SOURCE_DIR}/configuration.json
    COMMENT "Generating zsh completion _pheniqs"
    BYPRODUCTS "${ZSH_COMPLETION}"
)

# Add binary directory to include path for generated files
include_directories(${CMAKE_BINARY_DIR})

# Create the executable
add_executable(pheniqs 
    ${PHENIQS_SOURCES}
    "${VERSION_HEADER}"
    "${CONFIG_HEADER}"
)

# Link libraries and set properties
target_link_libraries(pheniqs PRIVATE ${PHENIQS_LIBRARIES} pthread)
target_include_directories(pheniqs PRIVATE ${HTSLIB_INCLUDE_DIR})

# Installation rules
install(TARGETS pheniqs DESTINATION ${BIN_PREFIX})

# Install zsh completion
install(FILES "${ZSH_COMPLETION}" 
        DESTINATION "${ZSH_PREFIX}/site-functions" 
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Add custom test targets
add_custom_target(test.api.configuration
    COMMAND ./test/api/configuration/run.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)

add_custom_target(test.api.illumina
    COMMAND ./test/api/illumina/run.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)

add_custom_target(test.api.io
    COMMAND ./test/api/io/run.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)

add_custom_target(test.api.prior
    COMMAND ./test/api/prior/run.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)

add_custom_target(test.pheniqs.BDGGG
    COMMAND ./test/BDGGG/run.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)

add_custom_target(test DEPENDS test.pheniqs.BDGGG)

# Clean targets for test results
add_custom_target(clean.test.api.configuration
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/test/api/configuration/result
)

add_custom_target(clean.test.api.illumina
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/test/api/illumina/result
)

add_custom_target(clean.test.api.io
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/test/api/io/result
)

add_custom_target(clean.test.api.prior
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/test/api/prior/result
)

add_custom_target(clean.test.pheniqs.BDGGG
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/test/BDGGG/result
)

# Help target
add_custom_target(pheniqs-help
    COMMAND ${CMAKE_COMMAND} -E echo "Building Pheniqs"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "CMake Targets:"
    COMMAND ${CMAKE_COMMAND} -E echo "  pheniqs              - Build the pheniqs executable (default)"
    COMMAND ${CMAKE_COMMAND} -E echo "  install              - Install pheniqs to CMAKE_INSTALL_PREFIX"
    COMMAND ${CMAKE_COMMAND} -E echo "  test                 - Run test suite"
    COMMAND ${CMAKE_COMMAND} -E echo "  clean                - Clean build artifacts"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Configuration Variables:"
    COMMAND ${CMAKE_COMMAND} -E echo "  CMAKE_INSTALL_PREFIX - Installation prefix (default: /usr/local)"
    COMMAND ${CMAKE_COMMAND} -E echo "  CMAKE_CXX_COMPILER   - C++ compiler to use"
    COMMAND ${CMAKE_COMMAND} -E echo "  BUILD_STATIC         - Build with static libraries (ON/OFF)"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Pheniqs depends on:"
    COMMAND ${CMAKE_COMMAND} -E echo "  zlib       - https://zlib.net"
    COMMAND ${CMAKE_COMMAND} -E echo "  bz2        - http://www.bzip.org"
    COMMAND ${CMAKE_COMMAND} -E echo "  xz/lzma    - https://tukaani.org/xz"
    COMMAND ${CMAKE_COMMAND} -E echo "  libdeflate - https://github.com/ebiggers/libdeflate (optional)"
    COMMAND ${CMAKE_COMMAND} -E echo "  htslib     - http://www.htslib.org"
    COMMAND ${CMAKE_COMMAND} -E echo "  rapidjson  - http://rapidjson.org"
)

# Print configuration
message(STATUS "")
message(STATUS "Pheniqs Configuration:")
message(STATUS "  Version:              ${PHENIQS_VERSION}")
message(STATUS "  Compiler:             ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ Standard:         ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install Prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Static Build:         ${BUILD_STATIC}")
message(STATUS "")
