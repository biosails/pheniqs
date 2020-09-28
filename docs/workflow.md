---
layout: default
title: "Workflow"
permalink: /workflow
id: workflow
---

* placeholder
{:toc}

Examples of how to configure Pheniqs for a handful of published experimental designs are provided below to help users get started with their own applications.

Different types of experiments use different barcoding schemes. The number, location, and type of barcodes to be extracted, which will vary depending on the specific experimental design. Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experiments. Different types of sequences are extracted using _tokens_ and handled using _transform patterns_.

See the [**workflow**](workflow) page for more information on common types of experimental designs, read layouts, tokenization, and transform patterns.

The [**Tokenization and Transform Patterns**](tokenization) page provides a high-level view of how transform patterns are constructed and the output they provide.

To configure Pheniqs for any particular workflow, it helps to have some idea about next-gen sequencing technology, adapters and barcode sets, Illumina's bcl2fastq software, FASTQ and SAM file formats, and some general terminology.


# Vignettes

Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. Follow the links below for descriptions of barcoding schemes for different experimental scenarios:

[Prior estimation with the Illumina Python API](vignettes/illumina_python_api)
: This vignette will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the Pheniqs python API. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors, and demultiplex the run. It loosely applies to almost every standard sample multiplex Illumina run.

[Standard Illumina sample demultiplexing](vignettes/illumina.html)
: This vignette will walk you through writing configuration files for demultiplexing a standard Illumina high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary#phred_adjusted_maximum_likelihood_decoding).

[Fluidigm with a sample and a cellular tag](vignettes/fluidigm.html)
: This vignette will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary#phred_adjusted_maximum_likelihood_decoding).

[Droplet Sequencing](vignettes/droplet.html)
: A popular way to perform single-cell RNA-seq profiling is to use microfluidics to generate individual droplets that each contain an individual cell and a bead carrying a cellular barcode. This page provides barcoding scenarios for several different droplet sequencing methods.

[Split-Pool Indexing](vignettes/split_pool.html)
: Another way to perform single-cell profiling is to perform a series of sequential cell pooling and splitting onto multi-well plates that are barcoded by row and column. Two or more rounds of barcoding in this manner will produce a unique cellular barcode for each cell that is composed of multiple independent elements. This page provides barcoding scenarios for published split-pool methods.

[Multimodal Profiling](vignettes/multimodal.html)
: Recently, methods have been developed to profile multiple layers of cellular information in the same experiment, such as mRNA abundance and protein levels for select cell-surface markers. This page provides barcoding scenarios for some of the multimodal profiling methods published so far.
