# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is Pheniqs

Pheniqs (PHilology ENcoder wIth Quality Statistics) is a multithreaded C++11 barcode classifier for high-throughput next-generation sequencing. It demultiplexes samples using probabilistic Bayesian classification and supports FASTQ, SAM/BAM/CRAM formats.

## Build Commands

### CMake (preferred, current branch is `cmake_conversion`)

```bash
mkdir build && cd build
cmake ..
make -j4
```

**Key CMake options:**
- `-DBUILD_STATIC=ON` — static binary
- `-DVENDOR_HTSLIB=ON` — build htslib locally (also set for its deps if vendoring any)
- `-DVENDOR_ZLIB=ON`, `-DVENDOR_BZIP2=ON`, `-DVENDOR_LZMA=ON`, `-DVENDOR_LIBDEFLATE=ON` — vendor compression libs
- `-DVENDOR_RAPIDJSON=ON` — vendor header-only RapidJSON
- `-DCMAKE_INSTALL_PREFIX=/path` — install prefix (default: `/usr/local`)

**Constraint:** When using `BUILD_STATIC=ON`, if any of htslib's deps are vendored, `VENDOR_HTSLIB` must also be `ON`.

Vendored dependencies are downloaded from the internet at build time and installed to `build/external/install/`.

### Makefile (legacy)

```bash
make -j4
make with-static=1 -j4   # static build
```

Dependencies must be installed at `$(PREFIX)` (default `/usr/local`).

### Generated files

Two headers are generated at build time into the build directory:
- `version.h` — from git describe via `cmake/GenerateVersionHeader.cmake`
- `configuration.h` — from `configuration.json` via `tool/pheniqs-configuration-api.py header`
- `_pheniqs` — zsh completion, also from `configuration.json`

## Running Tests

### Unit and integration tests (C++)

All tests live in `test/unit/` and run via a single binary `pheniqs_unit_tests`. Must be run from the repository root.

```bash
# from build dir
cmake --build build --target test     # build deps + run all tests
cmake --build build --target test.unit  # same (alias)
```

Test files:
- `test_phred.cpp` — `PhredScale` lookup table correctness
- `test_sequence.cpp` — `Sequence`, `ObservedSequence`, `Observation`
- `test_demux.cpp` — full demultiplexing pipeline (calls `Pipeline` in-process with BDGGG test data, checks SAM output at `test/BDGGG/result/demux_test_output.sam`)

**Known failing test:** `test_append_corrected_nonempty_dest` exposes a bug in `ObservedSequence::append_corrected` (`sequence.h:389`): when the destination already has content (`length > 0`), quality assignment uses `original.code[dest_length + i]` instead of `original.code[start + i]`, corrupting quality scores for uncorrected bases.

## Architecture

### Entry point and pipeline

`pheniqs.cpp` → `Pipeline` (pipeline.h/cpp) → reads config, manages job queue → `Job` (job.h/cpp) → `Transcode` (transcode.h/cpp) orchestrates reading, classifying, and writing reads.

### Classification algorithms (three decoders)

All decoders inherit from `Decoder<T>` template:

| File | Class | Algorithm |
|------|-------|-----------|
| `naive.h/cpp` | `NaiveMolecularDecoder` | Direct extraction, no classification |
| `mdd.h/cpp` | `MdDecoder<T>` | Minimum Distance (edit distance + quality masking) |
| `pamld.h/cpp` | `PamlDecoder<T>` | Probabilistic Accurate ML — Bayesian with confidence scores |

`PamlDecoder` is the primary, publication-described algorithm. `MdDecoder` is the fallback edit-distance approach.

### I/O layer

- `fastq.h/cpp` — FASTQ reader
- `hts.h/cpp` — HTSLib wrapper (BAM/SAM/CRAM)
- `feed.h/cpp` — abstract feed interface
- `proxy.h/cpp` — data source proxy

### Core data model

- `read.h/cpp` — read object
- `sequence.h/cpp` — DNA sequence + quality scores
- `barcode.h/cpp` — barcode definitions with concentration priors
- `auxiliary.h/cpp` — SAM auxiliary tags (RG, BC, QT, etc.)
- `phred.h/cpp` — Phred quality scale math and substitution lookup tables

### Configuration and CLI

- `interface.h/cpp` — CLI parsing
- `configuration.json` — canonical definition of all CLI options; used to generate `configuration.h` and zsh completion
- `json.h/cpp` — RapidJSON wrapper
- `atom.h/cpp` — JSON element abstraction

### Selector/transform pipeline

`selector.h/cpp` and `transform.h/cpp` implement a rule engine for extracting barcode segments from arbitrary positions in reads. `multiplex.h/cpp` handles the multiplexing/demultiplexing logic on top of decoded barcodes.

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| htslib | 1.23 | BAM/CRAM support (primary dep) |
| zlib | 1.3.1 | Compression |
| bzip2 | 1.0.8 | Compression |
| xz/liblzma | 5.8.2 | Compression |
| libdeflate | 1.15 | Optional, faster gzip |
| RapidJSON | 1.1.0 | Header-only JSON |
| pthread | system | Multithreading |
