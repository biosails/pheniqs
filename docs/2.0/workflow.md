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
        <li><a                  href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a class="active"   href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Build</a></li>
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

Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experimental designs. The conceptual framework of sequence classification and demultiplexing employed by Pheniqs is summarized below. Examples of how to configure Pheniqs for a handful of published experimental designs may be found in the [vignettes section](workflow.html) of the documentation.

# Experimental designs

Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. Some common designs for the Illumina platform are illustrated here:

![experimental designs](/pheniqs/assets/img/diagram8.png)

<a name="illumina_python_api" />
[Prior estimated Illumina with the python API](illumina_python_api.html)
: This vignette will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the Pheniqs python API. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors and demultiplex the run. It loosely applies to almost every standard sample multiplex Illumina run.

<a name="standard_illumina" />
[Standard Illumina sample demultiplexing](illumina.html)
: This vignette will walk you through writing configuration files for demultiplexing a standard Illumina high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
[Fluidigm with a sample and a cellular tag](fluidigm.html)
: This vignette will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

# Read anatomy

Illumina sequencing platforms typically produce four different sequence elements: two Index sequences, referred by Illumina as the **i5** and **i7** barcodes, and two Insert sequences, referred by Illumina as **read 1** and **read 2**. Collectively, these are referred to as read segments. For example, consider a [standard paired-end, dual index library design](illumina.html):

<!-- ![read anatomy](/pheniqs/assets/img/diagram1.png) -->

<img src="/pheniqs/assets/img/diagram1.png" class="figure_medium" />

<!-- <img src="/pheniqs/assets/img/diagram1.png" width="75%" /> -->

The read segments for this standard design thus comprise two technical sequences (referred by Illumina as I1, I2) and two biological sequences (referred by Illumina as R1, R2):

![read anatomy](/pheniqs/assets/img/diagram2.png)

# Sequence Classification

The combination of the barcodes contained in the **I1** and **I2** index positions specifies the sample library. With standard dual indexing, up to 96 distinct sample libraries can be pooled and run together in a single sequencing lane.

To identify which biological sequences belong to which library, the sequences belonging each one need to be separated from each other. This process of deconvolving libraries is called demultiplexing and is done by classifying each of the sequences using the barcode indexes:

<center>
![read anatomy](/pheniqs/assets/img/diagram5.png){ width=60% }
</center>

>**Note** The i5 adaptor sequences specified in sample sheets will be reverse complemented for platforms that read the I2 index on the bottom strand. The i7 sequences in sample sheets are always reverse complemented relative to the original adaptor sequences since they are read from the top strand
{: .example}

# Tokenization

Pheniqs uses [tokens](manual.html#tokenization) to reference and extract information from different read segments by specifying where to look for different classes of sequence elements (i.e. barcodes, biological sequences). Each element of interest is defined by an offset relative to the beginning of a given read segment (in this example I1, I2, R1, R2) and a length. It is important to note that Pheniqs uses [zero based](glossary.html#zero_based_coordinate) indexing, so the first read to come off the machine will be referred to as Segment 0, and so on:

<center>
![read anatomy](/pheniqs/assets/img/diagram7.png){ width=60% }
</center>

For a standard paired-end, dual indexed Illumina run, the sample barcodes usually comprise the full I1 and I2 read segments. Because Illumina sequencing is calibrated in relation to the previously sequenced base, those segments are sometimes sequenced one nucleotide longer than necessary to ensure good quality on the last nucleotide. The biological sequences start at the first position of R1 and R2 and extend for the full number of cycles run (typically 75, 100, or 150 nucleotides).

For this design, the barcode tokens begin at position 0 in I1 and I2 and extend for 10 bases. The tokens for biological sequences begin at position 0 of Read 1 and Read 2 and extend for the full span of those read segments.

# Input / Output

Pheniqs can manipulate [SAM, BAM and CRAM](glossary.html#htslib) files as well as uncompressed and gzip compressed [FASTQ](glossary.html#fastq). Configuration and reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration.

All of these can be specified with standard templates and may be overridden by additional directives within the configuration file.

After sequence classification, the various tokenized sequences extracted are written along with their respective confidence scores to specific [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf) within the [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format, as outlined below.

| Name                                      | Description                                                    | Example                       |
| :---------------------------------------- | :------------------------------------------------------------- | :---------------------------- |
| **[RG](glossary.html#rg_auxiliary_tag)**  | Read group identifier matching an RG filed in the header.      | H7LT2DSXX:1:GAACTGAGCGTCGTGGAGCG  |
| **[BC](glossary.html#bc_auxiliary_tag)**  | Raw uncorrected sample barcode sequence.                       | GAACTGAGCG-TCGTGGAGCG             |
| **[QT](glossary.html#qt_auxiliary_tag)**  | Phred quality of the sample barcode sequence in the BC tag.    | ,FF::F:F:F-,,FF::FF:F             |
| **[XB](glossary.html#xb_auxiliary_tag)**  | The probability that sample barcode decoding is incorrect.     | 2.27479e-06                   |
| **[CB](glossary.html#cb_auxiliary_tag)**  | Cellular identifier.                                           | ACTGCATA                      |
| **[CR](glossary.html#cr_auxiliary_tag)**  | Raw uncorrected cellular barcode sequence.                     | ACTGCATT                      |
| **[CY](glossary.html#cr_auxiliary_tag)**  | Phred quality of the cellular barcode sequence in the CR tag.  | ,,FF::FF                      |
| **[XC](glossary.html#xc_auxiliary_tag)**  | The probability that Cellular barcode decoding is incorrect.   | 2.27479e-06                   |
| **[MI](glossary.html#mi_auxiliary_tag)**  | Molecular Identifier.                                          |                               |
| **[RX](glossary.html#rx_auxiliary_tag)**  | Molecular barcode sequence, either corrected or uncorrected.   |                               |
| **[QX](glossary.html#qx_auxiliary_tag)**  | Phred quality of the molecular barcode sequence in the RX tag. |                               |
| **[OX](glossary.html#ox_auxiliary_tag)**  | Raw uncorrected molecular barcode sequence.                    |                               |
| **[BZ](glossary.html#bz_auxiliary_tag)**  | Phred quality of the molecular barcode sequence in the OX tag. |                               |
| **[XM](glossary.html#xm_auxiliary_tag)**  | The probability that molecular barcode decoding is incorrect.  |                               |
