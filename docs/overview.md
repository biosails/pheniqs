---
layout: default
title: "Overview"
permalink: /overview/
id: overview
---

* placeholder
{:toc}

The framework of a Pheniqs run includes three main steps: validating run parameters, extracting sequence elements from input sequence reads, and decoding each element according to a set of rules provided in the configuration:

![overview](/pheniqs/assets/img/pheniqs_overview_web2.png)

The essential components of a Pheniqs workflow are:
+ **Input**: FASTQ or SAM-formatted sequence files.
+ **Configuration**: All runtime directives, including I/O, sequence elements of interest, barcode sets, metadata, and any prior information about sample distributions. Pheniqs can be run with default parameters, but in most cases a configuration file will be required.
+ **Tokenization**: Each read segment is parsed to extract each sequence element to be processed, as specified in the configuration.
+ **Decoding**: For each barcode, Pheniqs will perform either probabilistic decoding (PAMLD, preferred) or simple minimum distance decoding (MDD), as specified in the configuration.
+ **Output**: Biological sequences, observed and inferred barcode sequences, quality scores, and decoding error probabilities are emitted as output. Sequence Alignment/Map (SAM) format is preferred, but FASTQ may also be emitted.
+ **Run Report**: Summary statistics about the decoding run are computed and written in a machine-readable JSON format, which can be easily parsed for visual display.

# Input

Pheniqs is designed to take **sequence** files as input. All standard sequence file formats are accepted. Most commonly, three or four FASTQ files, as emitted by Illumina sequencers, will be used as input. Pheniqs can manipulate both uncompressed and gzip compressed [FASTQ](glossary.html#fastq), as well as [SAM, BAM and CRAM](glossary.html#htslib) files.

# Configuration

Pheniqs also needs **configuration** directives in order to execute a decoding run. Required parameters include input / output file names and paths, barcode sets, transform directives for extracting tokens, and a variety of metadata, such as expected proportions of multiplexed sample libraries.

If no user-provided configuration file is available, Pheniqs will run with defaults for all parameters, but in most cases a configuration file will be needed to specify all of the run-specific information. Configuration files are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration with automated pipelines. A summary of all configuration directives is provided in the [Configuration](configuration.html) page.

At the beginning of a run, Pheniqs will **compile and validate** the configured parameters and will abort with explicit error messages upon any validation failure.

If **prior estimation** of barcode distributions are planned, a preliminary run will be executed and the configuration will be updated accordingly.

# Tokenization

Due to its flexible syntax for parsing read segments, Pheniqs can accommodate virtually any experimental design. Configuration directives will handle any combination of biological and technical sequences, such as barcoded multiplexed sample libraries, cellular indexes, and UMIs, for both bulk and single-cell experimental designs.

An overview of how Pheniqs parses sequence reads is provided in the [Tokens](tokenization.html) page.

Examples of how to configure tokenization transform patterns for a handful of published experimental designs may be found in the [vignettes section](vignettes.html) of the documentation.

# Decoding

Pheniqs currently implements two types of decoders to infer barcode sequences to be used for sequence classification:

+ A standard **minimum distance decoder (MDD)** that uses Hamming (edit) distance and simple string matching to allow zero or more errors per barcode, and
+ A **Phred-adjusted maximum likelihood decoder (PAMLD)**, which consults sequence quality scores and prior sample distributions to compute the full posterior probability for observed barcodes. PAMLD implements two successive filters to determine decoding success or failure, a _noise_ filter and a _high confidence_ filter:

![PAMLD](/pheniqs/assets/img/pamld.png)
{: .img-small}
<!-- <img src="/pheniqs/assets/img/pamld.png" style="img-small" /> -->

Reads with a lower conditional probability than random sequences fail the noise filter and are classified as noise without further consideration. Reads with a posterior probability that does not meet the confidence threshold fail the confidence filter; these reads are classified, but they are marked as "qc fail" so the confidence threshold can be reconsidered at alater stage. A full description of the mathematics behind Pheniqs, as well as performance evaluations and comparisons with other decoding methods, may be found [here]().

Pheniqs is designed to accommodate the addition of alternative decoders, which can be added as derived classes of a generic decoder object.

# Output

Biological sequence read segments are emitted along with observed and corrected barcode sequences and decoding error probability scores for each barcode. All of these may be specified with standard templates and may be overridden by additional directives within the configuration file.

The error score is computed as one minus the confidence score; for compound barcodes, it is one minus the product of the individual confidence scores. An example of the output is provided on the [Tokenization](tokenization.html) page.

# Run Report

Summary statistics for each run are also generated. Run reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded, which makes them easy to parse and use for producing tabular data and visualizations.
