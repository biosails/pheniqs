---
layout: default
title: "Pheniqs"
permalink: /
id: home
---

![transform patterns](/pheniqs/assets/img/transform_patterns.png)

Pheniqs is a flexible generic barcode classifier for high-throughput next-gen sequencing that caters to a wide variety of experimental designs and has been designed for efficient data processing. Features include:

### Powerful and intuitive syntax
- Addresses index tags in arbitrary locations along reads
- Easily handles any number of combinatorial barcode tags
- Classifies standard barcode types: Sample, Cellular, and Molecular Index
- Easily accommodates custom barcode types
- Directly writes barcodes to standard or custom BAM fields
- Eliminates need for pre- or post-processing of barcode tags

### Noise and quality aware probabilistic classifier
- Increases accuracy over standard edit distance methods
- Reports classification error probabilities in standard SAM auxiliary tags
- Modular design allows addition of new classifiers

### Robust engineering
- Multithreaded C++ implementation optimized for speed
- POSIX standard stream integration
- Directly interfaces with low level HTSLib C API
- Performance scales linearly with the number of available processing cores

### Easy to install or build
- Built-in package manager can build dependencies and binaries from scratch
- Stable releases available from Bioconda
- Available in a Docker container
- Portable compiled binaries available
- Easily installed on clusters or cloud without elevated permissions

### Easy to use
- Simple command line syntax with autocomplete
- Reusable configuration templates
- Preconfigured barcode library sets
- Reads and writes multiple file formats: FASTQ, SAM/BAM/CRAM
- Fast standalone file format interconversion
- Helper scripts assist in configuration file bootstraping
- Facilitates more robust and reproducible downstream analysis

On this website, you will find everything you need to get started with Pheniqs: Installation and configuration, classifying reads, vignettes that will walk you through processing popular experimental designs, and information about how to leverage standardized SAM auxiliary tags to improve the reproducibility of your published data.

Pheniqs runs on all modern POSIX systems and provides an easy to learn command line interface with autocomplete and an extensible reusable configuration syntax. Pheniqs is an ideal utility to pre- and post-process sequence reads for other bioinformatics tools, and it may also be used simply to rapidly and efficiently interconvert a variety of standard sequence file formats without invoking any of its barcode processing features.

For more advanced users and sequencing core managers, we provide detailed build instructions and a built-in package manager to easily build portable, [statically linked](glossary#static_linking), Pheniqs binaries for deployment on computing clusters. Developers can find code examples and API documention that enable them to expand Pheniqs with new classification algorithms and take advantage of the optimized multithreaded pipeline.

## Installing Pheniqs

You can install pheniqs using the [conda](install) package manager on Linux

>```shell
conda install -c bioconda -c conda-forge pheniqs
```

You can build the latest Pheniqs binary with this one line. Just paste it into a Linux or macOS terminal.

>```shell
bash -c "$(curl -fsSL https://raw.githubusercontent.com/biosails/pheniqs-build-api/master/install_pheniqs.sh)"
```
