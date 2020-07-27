---
layout: default
title: "Pheniqs"
permalink: /
id: home
---

Pheniqs is a generic DNA sequence multiplexer, demultiplexer and quality analyzer that caters to a wide variety of experimental designs. It is written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) and is can easily be to extend with new decoding algorithms. A python API for generating configuration files for common experimental designs is also available and actively maintained. Pheniqs is pronounced ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics.

## Pheniqs at a glance

**Format support:**
Pheniqs can manipulate [SAM, BAM and CRAM](/pheniqs/glossary#htslib) files as well as uncompressed and gzip compressed [FASTQ](/pheniqs/glossary#fastq). Configuration and reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration.

**Quality and prior aware DNA barcode decoding:**
in addition to the ubiquitous [minimum distance decoder](/pheniqs/glossary#minimum_distance_decoding) (**MDD**), Pheniqs offers a more accurate probabilistic [Phred-adjusted maximum likelihood decoder](/pheniqs/glossary#phred_adjusted_maximum_likelihood_decoding) (**PAMLD**). PAMLD consults base calling quality scores and a set of priors to compute the posterior probability of correctly decoding a barcode. The probability of an incorrect barcode assignment is reported in SAM auxiliary tags.

**Rich SAM metadata support:**
[JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files provide a powerful tool for manipulating the SAM [header](https://samtools.github.io/hts-specs/SAMv1.pdf) and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). During decoding Pheniqs populates the multiplex ([BC](/pheniqs/glossary#bc_auxiliary_tag) and [QT](/pheniqs/glossary#qt_auxiliary_tag)), molecular ([RX](/pheniqs/glossary#rx_auxiliary_tag), [QX](/pheniqs/glossary#qx_auxiliary_tag), [OX](/pheniqs/glossary#ox_auxiliary_tag), [BZ](/pheniqs/glossary#bz_auxiliary_tag) and [MI](/pheniqs/glossary#mi_auxiliary_tag)) and cellular (and [CB](/pheniqs/glossary#cb_auxiliary_tag)
[CR](/pheniqs/glossary#cr_auxiliary_tag) and [CY](/pheniqs/glossary#cr_auxiliary_tag)) standardized tags. Decoded sample identifiers are written to the standard [RG](/pheniqs/glossary#rg_auxiliary_tag) tag.

Barcode decoding error probabilities are encoded in tags reserved for local use until they are standardized: [XB](/pheniqs/glossary#xb_auxiliary_tag) for the multiplex, [XM](/pheniqs/glossary#xm_auxiliary_tag) for molecular and [XC](/pheniqs/glossary#xc_auxiliary_tag) for cellular. Pheniqs can also report the computed expected number of errors in the read in the [XE](/pheniqs/glossary#xe_auxiliary_tag) tag.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single, statically linked, compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and supports [zsh completion](https://en.wikipedia.org/wiki/Command-line_completion) for a more interactive command line experience.

**Easy to deploy:** Pheniqs depends on [HTSlib](http://www.htslib.org), [RapidJSON](http://rapidjson.org) and [zlib](https://zlib.net). HTSLib further depends on [bzip2](http://www.bzip.org), [LZMA](https://tukaani.org/xz) and optionally [libdeflate](https://github.com/ebiggers/libdeflate) for improved gzip compressed FASTQ manipulation. Pheniqs is regularly tested by [Travis](https://travis-ci.org/biosails/pheniqs) on Ubuntu and MacOS using most modern compilers. `pheniqs-build-api.py`, part of the Pheniqs python API, can build a portable, [statically linked](https://en.wikipedia.org/wiki/Static_library), pheniqs executable on most POSIX environments without elevated permissions by executing one simple command.

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0).

## Getting started
You can [install](https://biosails.github.io/pheniqs/2.0/install.html) Pheniqs from many common package managers or [build](https://biosails.github.io/pheniqs/2.0/build.html) it from source. The [quick tutorial](tutorial.md) will walk you though a very basic scenario of using Pheniqs. In the [manual](manual.md) you will find a more detailed discussion of each of the configuration parameters. In the [workflow](workflow.md) page you can find vignettes that will walk you through some standard protocols. **Please request specific protocols you wish to be included.**
