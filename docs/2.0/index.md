<!--
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
-->

<section id="navigation">
    <ul>
        <li><a class="active"   href="/pheniqs/">Home</a></li>
        <li><a                  href="/pheniqs/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Pheniqs
{:.page-title}

Pheniqs is a generic high throughput DNA sequence demultiplexer and quality analyzer written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11). Pheniqs is pronounced  ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics.

## Pheniqs at a glance

**File format support:** Pheniqs supports manipulating [SAM, BAM and CRAM](glossary.html#htslib) file formats as well as [FASTQ](glossary.html#fastq). Configuration and reports are encoded in [JSON](https://en.wikipedia.org/wiki/JSON) for easy portability and integration.

**Quality and prior aware barcode decoding:** In addition to the widespread [minimum distance decoder](glossary.html#minimum_distance_decoding) Pheniqs can decode barcodes using a [Phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that consults base calling quality scores and a set of user provided priors on the barcode prevalence distribution. The decoder reports accurately computed barcode decoding error probabilities in SAM auxiliary tags.

**Rich SAM metadata support:** Pheniqs configuration files provide a powerful tool for manipulating [SAM header](https://samtools.github.io/hts-specs/SAMv1.pdf) and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). During decoding Pheniqs writes multiplex ([BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag)), molecular ([RX](glossary.html#rx_auxiliary_tag), [QX](glossary.html#qx_auxiliary_tag), [OX](glossary.html#ox_auxiliary_tag), [BZ](glossary.html#bz_auxiliary_tag) and [MI](glossary.html#mi_auxiliary_tag)) and cellular (and [CB](glossary.html#cb_auxiliary_tag)
[CR](glossary.html#cr_auxiliary_tag) and [CY](glossary.html#cr_auxiliary_tag)) barcodes in the corresponding standardized tags. Decoded multiplex barcode identifiers are written to the standard [RG](glossary.html#rg_auxiliary_tag) tag, while the decoder declaration in the configuration file is used to further define read group header annotations. Barcode decoding error probabilities are encoded in tags reserved for local use until they are standardized: [XB](glossary.html#bc_auxiliary_tag) for the multiplex, [XM](glossary.html#xm_auxiliary_tag) for molecular and [XC](glossary.html#cr_auxiliary_tag) for cellular. Pheniqs can also report the computed expected number of errors in the read in the [XE](glossary.html#xe_auxiliary_tag) tag.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single, statically linked, compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and supports [zsh completion]({{ site.github.repository_url }}/blob/master/zsh/_pheniqs) for a more interactive command line experience.

**Easy to deploy:** Pheniqs depends on [HTSlib](http://www.htslib.org), [RapidJSON](http://rapidjson.org) and [zlib](https://zlib.net). HTSLib further depends on [bzip2](http://www.bzip.org), [LZMA](https://tukaani.org/xz) and optionally [libdeflate](https://github.com/ebiggers/libdeflate) for improved gzip compressed FASTQ manipulation. Pheniqs is regularly tested by [Travis](https://travis-ci.org/biosails/pheniqs) on Ubuntu and MacOS using most modern compilers. The provided pheniqs-tools utility can build a portable, [statically linked](https://en.wikipedia.org/wiki/Static_library), pheniqs executable on most POSIX environments without elevated permissions by executing one simple command.

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0).

## Getting started
To get you started, the [quick tutorial](tutorial.md) will walk you though a basic scenario of using Pheniqs. In the [manual](manual.md) you will find a more detailed discussion of each of the configuration parameters. Finally the [workflow](workflow.md) section demonstrates more involved read manipulation.
