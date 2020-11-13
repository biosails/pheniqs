# Pheniqs

For documentation check out the **[Pheniqs website](http://biosails.github.io/pheniqs)**.

Pheniqs is a generic DNA sequence multiplexer, demultiplexer and quality analyzer that caters to a wide variety of experimental designs. It is written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) and is can easily be to extend with new decoding algorithms. A python API for generating configuration files for common experimental designs is also available and actively maintained. Pheniqs is pronounced ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics.

## Pheniqs at a glance

**Format support:**
Pheniqs can manipulate [SAM, BAM and CRAM](https://biosails.github.io/pheniqs/2.0/glossary.html#htslib) files as well as uncompressed and gzip compressed [FASTQ](https://biosails.github.io/pheniqs/2.0/glossary.html#fastq). Configuration and reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration.

**Quality and prior aware DNA barcode decoding:**
in addition to the ubiquitous [minimum distance decoder](https://biosails.github.io/pheniqs/2.0/glossary.html#minimum_distance_decoding) (**MDD**), Pheniqs offers a more accurate probabilistic [Phred-adjusted maximum likelihood decoder](https://biosails.github.io/pheniqs/2.0/glossary.html#phred_adjusted_maximum_likelihood_decoding) (**PAMLD**). PAMLD consults base calling quality scores and a set of priors to compute the posterior probability of correctly decoding a barcode. The probability of an incorrect barcode assignment is reported in SAM auxiliary tags.

**Rich SAM metadata support:**
[JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files provide a powerful tool for manipulating the SAM [header](https://samtools.github.io/hts-specs/SAMv1.pdf) and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). During decoding Pheniqs populates the multiplex ([BC](https://biosails.github.io/pheniqs/2.0/glossary.html#bc_auxiliary_tag) and [QT](https://biosails.github.io/pheniqs/2.0/glossary.html#qt_auxiliary_tag)), molecular ([RX](https://biosails.github.io/pheniqs/2.0/glossary.html#rx_auxiliary_tag), [QX](https://biosails.github.io/pheniqs/2.0/glossary.html#qx_auxiliary_tag), [OX](https://biosails.github.io/pheniqs/2.0/glossary.html#ox_auxiliary_tag), [BZ](https://biosails.github.io/pheniqs/2.0/glossary.html#bz_auxiliary_tag) and [MI](https://biosails.github.io/pheniqs/2.0/glossary.html#mi_auxiliary_tag)) and cellular (and [CB](https://biosails.github.io/pheniqs/2.0/glossary.html#cb_auxiliary_tag)
[CR](https://biosails.github.io/pheniqs/2.0/glossary.html#cr_auxiliary_tag) and [CY](https://biosails.github.io/pheniqs/2.0/glossary.html#cr_auxiliary_tag)) standardized tags. Decoded sample identifiers are written to the standard [RG](https://biosails.github.io/pheniqs/2.0/glossary.html#rg_auxiliary_tag) tag.

Barcode decoding error probabilities are encoded in tags reserved for local use until they are standardized: [XB](https://biosails.github.io/pheniqs/2.0/glossary.html#xb_auxiliary_tag) for the multiplex, [XM](https://biosails.github.io/pheniqs/2.0/glossary.html#xm_auxiliary_tag) for molecular and [XC](https://biosails.github.io/pheniqs/2.0/glossary.html#xc_auxiliary_tag) for cellular. Pheniqs can also report the computed expected number of errors in the read in the [XE](https://biosails.github.io/pheniqs/2.0/glossary.html#xe_auxiliary_tag) tag.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single, statically linked, compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and supports [zsh completion](https://en.wikipedia.org/wiki/Command-line_completion) for a more interactive command line experience.

**Easy to deploy:** Pheniqs depends on [HTSlib](http://www.htslib.org), [RapidJSON](http://rapidjson.org) and [zlib](https://zlib.net). HTSLib further depends on [bzip2](http://www.bzip.org), [LZMA](https://tukaani.org/xz) and optionally [libdeflate](https://github.com/ebiggers/libdeflate) for improved gzip compressed FASTQ manipulation. Pheniqs is regularly tested by [Travis](https://travis-ci.org/biosails/pheniqs) on Ubuntu and MacOS using most modern compilers. `pheniqs-build-api.py`, part of the Pheniqs python API, can build a portable, [statically linked](https://en.wikipedia.org/wiki/Static_library), pheniqs executable on most POSIX environments without elevated permissions by executing one simple command.

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0).

## Getting started
You can [install](https://biosails.github.io/pheniqs/2.0/install.html) Pheniqs from many common package managers or [build](https://biosails.github.io/pheniqs/2.0/build.html) it from source. The [quick tutorial](https://biosails.github.io/pheniqs/2.0/tutorial.html) will walk you though a very basic scenario of using Pheniqs. In the [manual](https://biosails.github.io/pheniqs/2.0/manual.html) you will find a more detailed discussion of each of the configuration parameters. In the [workflow](https://biosails.github.io/pheniqs/2.0/workflow.html) page you can find vignettes that will walk you through some standard protocols. **Please request specific protocols you wish to be included.**

## Getting started
Head over to the [Pheniqs web site](https://biosails.github.io/pheniqs) for more details.

[![Build Status](https://travis-ci.org/biosails/pheniqs.svg?branch=master)](https://travis-ci.org/biosails/pheniqs)
