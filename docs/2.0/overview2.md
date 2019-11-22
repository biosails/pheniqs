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
        <li><a                  href="/pheniqs/2.0/">Home</a></li>
        <li><a                  href="/pheniqs/2.0/quickstart.html">Quickstart</a></li>
        <li><a class="active"   href="/pheniqs/2.0/workflow.html">Overview</a></li>
        <li><a                  href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a class="active"   href="/pheniqs/2.0/transform.html">Tokens</a></li>
        <li><a class="active"   href="/pheniqs/2.0/vignettes.html">Vignettes</a></li>
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Overview
{:.page-title}

* placeholder
{:toc}

The essential components of a Pheniqs workflow are illustrated below:

![overview](/pheniqs/assets/img/Pheniqs_overview_web.png)

+ Input files: FASTQ or SAM-formatted files.
+ Configuration file: Pheniqs can be run with default parameters, but in most cases a configuration file will be required.
+ Tokenization: Each read segment is parsed to extract the sequences of interest as specified in the configuration.
+ Decoding: For each barcode, Pheniqs may perform either probabilistic decoding (PAMLD, preferred) or simple minimum distance decoding (MDD).
+ Output: Biological sequences, observed and inferred barcode sequences, quality scores, and decoding error probabilities are emitted as output. Sequence Alignment/Map (SAM) format is preferred, but FASTQ may also be emitted.
+ Run Report: Summary statistics about the decoding run are also provided.

## Input

Pheniqs is designed to take sequence files and configuration directives as input. All standard formats are accepted. Most commonly, three or four FASTQ files, as emitted by Illumina sequencers, will be used as input.

Pheniqs can manipulate [SAM, BAM and CRAM](glossary.html#htslib) files as well as uncompressed and gzip compressed [FASTQ](glossary.html#fastq).

## Configuration

Pheniqs needs to know a variety of things before it can proceed with a decoding run. These include file paths and names of input / output files, barcode sets, transform directives for extracting tokens, and a variety of metadata, such as expected sample proportions of different libraries that have been multiplexed.

Configuration files are are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration with automated pipelines. A summary of all configuration directives is provided in the [Configuration](configuration.html) page.

## Tokenization

Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experimental designs.

Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. An overview of how Pheniqs parses sequence reads is provided in the [Tokens](transform.html) page.

Examples of how to configure tokenization transform patterns for a handful of published experimental designs may be found in the [vignettes section](vignettes.html) of the documentation.

## Decoding

## Output

... see Tokenization page ...

All of these can be specified with standard templates and may be overridden by additional directives within the configuration file.

## Run Report

... Reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded, which makes them easy to parse and use for producing tabular data and visualizations.
