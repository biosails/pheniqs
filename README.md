# Pheniqs

**NOTICE: Documentation for pheniqs 2.0 is undergoing and is expected to be avialbale in the next few days.**

Documentation and examples can be found on the **[Pheniqs website](http://biosails.github.io/pheniqs)**.

Pheniqs is a generic high throughput DNA sequence multiplexer/demultiplexer and quality analyzer written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11). Pheniqs is pronounced  ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics.

## Pheniqs at a glance

**Format support:**
Pheniqs can manipulate [SAM, BAM and CRAM](glossary.html#htslib) files as well as uncompressed and gzip compressed [FASTQ](glossary.html#fastq). Configuration and reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded.

**Quality and prior aware DNA barcode decoding:**
Pheniqs offers a choice of barcode decoding algorithms: The naive [minimum distance decoder](glossary.html#minimum_distance_decoding) (**MDD**) and a more advanced probabilistic [Phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) (**PAMLD**). PAMLD consults base calling quality scores and a set of user provided priors on the barcode prevalence distribution to compute an accurate decoding probability. The probability of an incorrect barcode assignment is reported in SAM auxiliary tags.

**Rich SAM metadata support:**
Pheniqs configuration files provide a powerful tool for manipulating the SAM [header](https://samtools.github.io/hts-specs/SAMv1.pdf) and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). During decoding Pheniqs populates the multiplex ([BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag)), molecular ([RX](glossary.html#rx_auxiliary_tag), [QX](glossary.html#qx_auxiliary_tag), [OX](glossary.html#ox_auxiliary_tag), [BZ](glossary.html#bz_auxiliary_tag) and [MI](glossary.html#mi_auxiliary_tag)) and cellular (and [CB](glossary.html#cb_auxiliary_tag)
[CR](glossary.html#cr_auxiliary_tag) and [CY](glossary.html#cr_auxiliary_tag)) standardized tags. Decoded multiplex barcode identifiers are written to the standard [RG](glossary.html#rg_auxiliary_tag) tag. The decoder declaration in the configuration file can be used to further manipulate read group header annotations. Barcode decoding error probabilities are encoded in tags reserved for local use until they are standardized: [XB](glossary.html#bc_auxiliary_tag) for the multiplex, [XM](glossary.html#xm_auxiliary_tag) for molecular and [XC](glossary.html#cr_auxiliary_tag) for cellular. Pheniqs can also report the computed expected number of errors in the read in the [XE](glossary.html#xe_auxiliary_tag) tag.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single, statically linked, compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and supports [zsh completion]({{ site.github.repository_url }}/blob/master/zsh/_pheniqs) for a more interactive command line experience.

**Easy to deploy:** Pheniqs depends on [HTSlib](http://www.htslib.org), [RapidJSON](http://rapidjson.org) and [zlib](https://zlib.net). HTSLib further depends on [bzip2](http://www.bzip.org), [LZMA](https://tukaani.org/xz) and optionally [libdeflate](https://github.com/ebiggers/libdeflate) for improved gzip compressed FASTQ manipulation. Pheniqs is regularly tested by [Travis](https://travis-ci.org/biosails/pheniqs) on Ubuntu and MacOS using most modern compilers. The provided pheniqs-tools utility can build a portable, [statically linked](https://en.wikipedia.org/wiki/Static_library), pheniqs executable on most POSIX environments without elevated permissions by executing one simple command.

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0).

## Getting started

Head over to the [Pheniqs web site](http://biosails.github.io/pheniqs) for a [quick tutorial](http://biosails.github.io/pheniqs/tutorial.html), [documentation](http://biosails.github.io/pheniqs/manual.html) and [building instructions](http://biosails.github.io/pheniqs/building.html).

[![Build Status](https://travis-ci.org/biosails/pheniqs.svg?branch=master)](https://travis-ci.org/biosails/pheniqs)
