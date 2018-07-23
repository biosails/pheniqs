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
        <li><a class="active"   href="/pheniqs/1.0/">Home</a></li>
        <li><a                  href="/pheniqs/1.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/1.0/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/1.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/1.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/1.0/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/1.0/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Pheniqs
{:.page-title}

Pheniqs is a generic high throughput DNA sequence demultiplexer and quality analyzer written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11). Pheniqs is pronounced  ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics. 

## Pheniqs at a glance

**Flexible file formats:** Pheniqs supports [FASTQ](glossary.html#fastq) as well as the [SAM, BAM and CRAM](glossary.html#htslib) file formats.

**Quality aware barcode decoding:** In addition to the widespread [minimum distance decoder](glossary.html#minimum_distance_decoding) Pheniqs introduces the [Phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that consults base calling quality scores and estimate the barcode decoding error probability for each read.

**Standardized tags:** Pheniqs encodes the raw multiplex barcode sequence and quality in the corresponding [BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag) [SAM auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf), as well as classifying multiplexed read with the [RG](glossary.html#rg_auxiliary_tag) auxiliary tag.

**Community tags:** Pheniqs also encodes the raw molecular barcode sequence and quality in the corresponding community adopted [RX](glossary.html#rx_auxiliary_tag) and [QX](glossary.html#qx_auxiliary_tag) auxiliary tags.

**Proposed tags:** Pheniqs proposes to standardize three new auxiliary tags: [DQ](glossary.html#dq_auxiliary_tag), [PX](glossary.html#px_auxiliary_tag) and [EE](glossary.html#ee_auxiliary_tag) for encoding the multiplex barcode decoding error probability, molecular barcode decoding error probability and the expected number of errors, respectively.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and even comes with its own [zsh completion]({{ site.github.repository_url }}/blob/master/zsh/_pheniqs) script for a more interactive command line experience. 

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Few dependencies:** Pheniqs depends only on [HTSlib](https://github.com/samtools/htslib) and [RapidJSON](https://github.com/miloyip/rapidjson). Both are widely packaged, easy to install from source and maintained by a highly active community.

**Open Source:** Pheniqs is released under the terms of the [AGPL 3.0 license agreement](http://opensource.org/licenses/AGPL-3.0). 

## Getting started

To get you started, the [quick tutorial](tutorial.md) will walk you though a basic scenario of using Pheniqs. In the [manual](manual.md) you will find a more detailed discussion of each of the configuration parameters. Finally the [workflow](workflow.md) section demonstrates more involved read manipulation.

## Building

Pheniqs is developed on MacOS and Ubuntu and is tested to build with the supplied [Makefile]({{ site.github.repository_url }}/blob/master/Makefile) but should build on most POSIX systems. For more details head over to the [building](building.md) section.

