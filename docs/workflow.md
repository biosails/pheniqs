---
layout: default
title: "Workflow"
permalink: /workflow/
id: workflow
---

* placeholder
{:toc}

# Experimental designs

![read anatomy](/pheniqs/assets/img/diagram8.png)
<a name="illumina_python_api" />
[Prior estimated Illumina with the python API](illumina_python_api.html)
: This vignette will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the Pheniqs python API. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors and demultiplex the run. It loosely applies to almost every standard sample multiplex Illumina run.

<a name="standard_illumina" />
[Standard Illumina sample demultiplexing](illumina.html)
: This vignette will walk you through writing configuration files for demultiplexing a standard Illumina high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
[Fluidigm with a sample and a cellular tag](fluidigm.html)
: This vignette will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).
