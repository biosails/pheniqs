# Pheniqs

Pheniqs is a generic high throughput sequencing demultiplexer and quality analyzer written in multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11). Pheniqs is pronounced  ***phoe·nix*** and stands for **PH**ilology **EN**coder w**I**th **Q**uality **S**tatistics. 

Documentation and exmples are on the **[Pheniqs website](http://gunsaluspiano.github.io/pheniqs)**.

## Pheniqs at a glance

**Flexible file formats:** Pheniqs supports [FASTQ](http://gunsaluspiano.github.io/pheniqs/glossary.html#fastq) as well as the [SAM, BAM and CRAM](http://gunsaluspiano.github.io/pheniqs/glossary.html#htslib) file formats.

**Quality aware barcode decoding:** In addition to the widespread [minimum distance decoder](http://gunsaluspiano.github.io/pheniqs/glossary.html#minimum_distance_decoding) Pheniqs introduces the [Phred-adjusted maximum likelihood decoder](http://gunsaluspiano.github.io/pheniqs/glossary.html#phred_adjusted_maximum_likelihood_decoding) that consults base calling quality scores and estimate the barcode decoding error probability for each read.

**Standardized tags:** Pheniqs supports encoding the raw multiplex barcode sequence and quality in the corresponding [BC](http://gunsaluspiano.github.io/pheniqs/glossary.html#bc_auxiliary_tag) and [QT](http://gunsaluspiano.github.io/pheniqs/glossary.html#qt_auxiliary_tag) [SAM auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf), as well as classifying multiplexed read with the [RG](http://gunsaluspiano.github.io/pheniqs/glossary.html#rg_auxiliary_tag) auxiliary tag.

**Community tags:** Pheniqs also supports encoding the raw molecular barcode sequence and quality in the corresponding community adopted [RX](http://gunsaluspiano.github.io/pheniqs/glossary.html#rx_auxiliary_tag) and [QX](http://gunsaluspiano.github.io/pheniqs/glossary.html#qx_auxiliary_tag) auxiliary tags.

**Proposed tags:** Pheniqs proposes to standardize three new auxiliary tags: [DQ](http://gunsaluspiano.github.io/pheniqs/glossary.html#dq_auxiliary_tag), [PX](http://gunsaluspiano.github.io/pheniqs/glossary.html#px_auxiliary_tag) and [EE](http://gunsaluspiano.github.io/pheniqs/glossary.html#ee_auxiliary_tag) for encoding the multiplex barcode decoding error probability, molecular barcode decoding error probability and the expected number of errors, respectively.

**Command line friendly:** Written in multi threaded C++ Pheniqs can be deployed as a single compiled binary executable. Pheniqs follows the [POSIX command line syntax](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html) and even comes with its own [zsh completion](zsh_completion/_pheniqs) script for a more interactive command line experience. 

**Pipeline friendly:** [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration files and reports that can be streamed through [standard POSIX streams](https://en.wikipedia.org/wiki/Standard_streams) enable streamlined pipeline and [RESTful](https://en.wikipedia.org/wiki/Representational_state_transfer) service integration.

**Few dependencies:** Pheniqs depends only on [HTSlib](https://github.com/samtools/htslib) and [RapidJSON](https://github.com/miloyip/rapidjson). Both are widely packaged, easy to install from source and maintained by a highly active community.

**Open Source:** Pheniqs is released under the terms of the GNU General Public License. 

## Getting started
Head over to the [Pheniqs web site](http://gunsaluspiano.github.io/pheniqs) for a [quick tutorial](http://gunsaluspiano.github.io/pheniqs/tutorial.html), [documentation](http://gunsaluspiano.github.io/pheniqs/manual.html) and [building instructions](http://gunsaluspiano.github.io/pheniqs/building.html).
