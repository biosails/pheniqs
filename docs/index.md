---
layout: default
title: "Pheniqs"
permalink: /
id: home
---

Fast and accurate sequence manipulation

![transform patterns](/pheniqs/assets/img/transform_patterns.png)

Pheniqs is a generic high throughput next generation barcode classifier that caters to a wide variety of experimental designs. An intuitive addressing syntax allows researchers to easily classify reads by standard Illumina sample barcodes, multiple combinatorial cellular barcodes and molecular index tags in arbitrary locations along the read, all without pre or post processing their data. Pheniqs classifies reads using a noise and quality aware probabilistic classifier that is more accurate than widespread edit distance methods. Furthermore, by reporting the error probability for each classification in standasrd SAM auxiliary tags Pheniqs enables more robust and reproducible downstream analysis. To handle the rapid increase in sequencing throughput, a fine tuned multithreaded C++ implemenentation that directly interfaces with the low level HTSLib c api offers performance that scales linearly with the number of available processing cores.

Pheniqs runs on all modern POSIX systems, provides an easy to learn command line interface with autocomplete and extensible reusable configuration syntax. On this website, you will find everything you need to get started with Pheniqs: Installion and configuration, classifying reads, pre and post processing sequene data for other bioinformatics tools, vignettes that will walk you through processing popular experimental designs, and loads of information about how to leverage standardized SAM auxiliary tags to improve the reproduciblity of your published data.

For more advanced users and sequencing core managers we provide detailed build instructions and a built in packaghe manager that allows to easily build portable, statically linked, Pheniqs binaries for deployment on computing clusters. For developers we make available code examples and a documented API that allow to expand Pheniqs with new classification algorithms.

## installing Pheniqs

You install pheniqs using the conda package manager on Linux

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

Or build a portable binary from source with the built in package manager with
>```shell
pending shell script that can check out code and build Pheniqs with the built in package manager in one command 
```
