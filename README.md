# Pheniqs

Pheniqs is a generic high throughput DNA sequence demultiplexer and quality analyzer written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11). Pheniqs is pronounced  ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics.

Documentation and examples can be found on the **[Pheniqs website](http://biosails.github.io/pheniqs)**.

## Pheniqs at a glance

**Latest file formats:** Pheniqs supports [FASTQ](http://biosails.github.io/pheniqs/glossary.html#fastq) and the [Sequence Alignment/Map Formats](http://biosails.github.io/pheniqs/glossary.html#htslib) SAM, BAM and CRAM.

**Quality aware barcode decoding:** In addition to the widespread [minimum distance decoder](http://biosails.github.io/pheniqs/glossary.html#minimum_distance_decoding) Pheniqs introduces the more accurate [Phred-adjusted maximum likelihood decoder](http://biosails.github.io/pheniqs/glossary.html#phred_adjusted_maximum_likelihood_decoding) that consults base calling quality scores and estimate the barcode decoding error probability for each read.

**Standardized tags:** Pheniqs encodes the raw multiplex barcode sequence and quality in the [BC](http://biosails.github.io/pheniqs/glossary.html#bc_auxiliary_tag) and [QT](http://biosails.github.io/pheniqs/glossary.html#qt_auxiliary_tag) [SAM auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf), as well as classifying multiplexed read with the [RG](http://biosails.github.io/pheniqs/glossary.html#rg_auxiliary_tag) auxiliary tag.

**Community tags:** Pheniqs also encodes the raw molecular barcode sequence and quality in the corresponding community adopted [RX](http://biosails.github.io/pheniqs/glossary.html#rx_auxiliary_tag) and [QX](http://biosails.github.io/pheniqs/glossary.html#qx_auxiliary_tag) auxiliary tags.

**Proposed tags:** Pheniqs proposes to standardize three new auxiliary tags: [XB](http://biosails.github.io/pheniqs/glossary.html#dq_auxiliary_tag), [PX](http://biosails.github.io/pheniqs/glossary.html#px_auxiliary_tag) and [EE](http://biosails.github.io/pheniqs/glossary.html#ee_auxiliary_tag) for encoding the multiplex barcode decoding error probability, molecular barcode decoding error probability and the expected number of errors.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and even comes with its own [zsh completion](zsh_completion/_pheniqs) script for a more interactive command line experience.

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Few dependencies:** Pheniqs depends only on [HTSlib](https://github.com/samtools/htslib) and [RapidJSON](https://github.com/miloyip/rapidjson). Both are widely packaged, easy to install from source and maintained by a highly active community.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0).

## Getting started

Head over to the [Pheniqs web site](http://biosails.github.io/pheniqs) for a [quick tutorial](http://biosails.github.io/pheniqs/tutorial.html), [documentation](http://biosails.github.io/pheniqs/manual.html) and [building instructions](http://biosails.github.io/pheniqs/building.html).

[![Build Status](https://travis-ci.org/biosails/pheniqs.svg?branch=master)](https://travis-ci.org/biosails/pheniqs)
