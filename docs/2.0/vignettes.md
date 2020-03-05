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
        <li><a class="active"   href="/pheniqs/2.0/overview.html">Overview</a></li>
        <li><a                  href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Vignettes
{:.page-title}

* placeholder
{:toc}

Examples of how to configure Pheniqs for a handful of published experimental designs are provided below to help users get started with their own applications.

Different types of experiments use different barcoding schemes. The number, location, and type of barcodes to be extracted, which will vary depending on the specific experimental design. Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experiments. Different types of sequences are extracted using _tokens_ and handled using _transform patterns_.

See the [**Overview**](overview.html) page for more information on common types of experimental designs, read layouts, tokenization, and transform patterns.

The [**Tokenization and Transform Patterns**](transform.html) page provides a high-level view of how transform patterns are constructed and the output they provide.

To configure Pheniqs for any particular workflow, it helps to have some idea about next-gen sequencing technology, adapters and barcode sets, Illumina's bcl2fastq softare, FASTQ and SAM file formats, and some general terminology. Please see the [**Background**](background.html) page for some tips on these if you are not familiar with these.


# Vignettes

Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. Follow the links below for descriptions of barcoding schemes for different experimental scenarios:

<a name="illumina_python_api" />
![read anatomy](/pheniqs/assets/img/diagram8.png)
[Prior estimation with the Illumina Python API](vignettes/illumina_python_api.html)
: This vignette will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the Pheniqs python API. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors, and demultiplex the run. It loosely applies to almost every standard sample multiplex Illumina run.

<a name="standard_illumina" />
[Standard Illumina sample demultiplexing](vignettes/illumina.html)
: This vignette will walk you through writing configuration files for demultiplexing a standard Illumina high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
[Fluidigm with a sample and a cellular tag](vignettes/fluidigm.html)
: This vignette will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="droplet" />
[Droplet Sequencing](vignettes/droplet.html)
: A popular way to perform single-cell RNA-seq profiling is to use microfluidics to generate individual droplets that each contain an individual cell and a bead carrying a cellular barcode. This page provides barcoding scenarios for several different droplet sequencing methods.

<a name="split-pool" />
[Split-Pool Indexing](vignettes/split_pool.html)
: Another way to perform single-cell profiling is to perform a series of sequential cell pooling and splitting onto multi-well plates that are barcoded by row and column. Two or more rounds of barcoding in this manner will produce a unique cellular barcode for each cell that is composed of multiple independent elements. This page provides barcoding scenarios for published split-pool methods.

<a name="multimodal" />
[Multimodal Profiling](vignettes/multimodal.html)
: Recently, methods have been developed to profile multiple layers of cellular information in the same experiment, such as mRNA abundance and protein levels for select cell-surface markers. This page provides barcoding scenarios for some of the multimodal profiling methods published so far.