---
layout: default
title: "Pheniqs"
permalink: /
id: home
---

Fast and accurate sequence manipulation

![transform patterns](/pheniqs/assets/img/transform_patterns.png)

Pheniqs is a generic high throughput barcode classifier that caters to a wide variety of experimental designs. An intuitive addressing syntax allows researchers to easily classify reads by standard Illumina sample barcodes, multiple combinatorial cellular barcodes and molecular index tags in arbitrary locations along the read, all without pre or post processing their data. Pheniqs classifies reads using a noise and quality aware probabilistic classifier that is more accurate than widespread edit distance methods. Furthermore, by reporting the error probability for each classification in standard SAM auxiliary tags, Pheniqs enables more robust and reproducible downstream analysis. To handle the rapid increase in sequencing throughput, a fine tuned multithreaded C++ implementation, that directly interfaces with the low level HTSLib C API, offers performance that scales linearly with the number of available processing cores.

Pheniqs runs on all modern POSIX systems, provides an easy to learn command line interface with autocomplete and an extensible reusable configuration syntax. On this website, you will find everything you need to get started with Pheniqs: Installation and configuration, classifying reads, pre and post processing sequence reads for other bioinformatics tools, vignettes that will walk you through processing popular experimental designs, and information about how to leverage standardized SAM auxiliary tags to improve the reproducibility of your published data.

For more advanced users and sequencing core managers we provide detailed build instructions and a built in package manager that allows to easily build portable, statically linked, Pheniqs binaries for deployment on computing clusters. Developers can find code examples and API documention that enable them to expand Pheniqs with new classification algorithms and take advantage of the optimized multithreaded pipeline.

## installing Pheniqs

You can build the latest Pheniqs binary with this one line. Just paste it into a Linux or macOS terminal.

>```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/biosails/pheniqs-build-api/master/install_pheniqs.sh)"
```

You can also install pheniqs using the conda package manager on Linux

>```shell
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh
conda install -c bioconda -c conda-forge pheniqs
```

Or MacOS
>```shell
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-MacOSX-x86_64.sh
bash Anaconda3-5.0.1-MacOSX-x86_64.sh
conda install -c bioconda -c conda-forge pheniqs
```
