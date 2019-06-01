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
        <li><a                  href="/Pheniqs/2.0/">Home</a></li>
        <li><a                  href="/Pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/Pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/Pheniqs/2.0/build.html">Build</a></li>
        <li><a class="active"   href="/Pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/Pheniqs/2.0/best_practices.html">Best Practice</a></li>
        <li><a                  href="/Pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/Pheniqs/2.0/manual.html">Manual</a></li>
        <li><a                  href="/Pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/Pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Single index Fluidigm
{:.page-title}

This tutorial will walk you through demultiplexing a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The read has 3 segments, 1 biological from the DNA or RNA fragment and 2 technical containing the i7 sample index and a cell specific tag. If the results are written to SAM, BAM or CRAM the multiplex barcode and its quality scores are written to the [BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag) respectively, and the decoding error probabilities is written to the [XB](glossary.html#xb_auxiliary_tag).

## Input Read Layout

<div class="read" id="single_index">
  <div class="binding_primer p5">P5</div>
  <div class="sequencing_primer">SP5</div>
  <div class="cellular_barcode">tag</div>
  <div class="insert">insert</div>
  <div class="sequencing_primer">SP7</div>
  <div class="sample_barcode">i7</div>
  <div class="binding_primer p7">P7</div>
  <div class="clear"></div>

  <p><strong>P5</strong> and <strong>P7</strong> are the primers binding to the instrument. <strong>SP5</strong> and <strong>SP7</strong> are the sequencing primers. <strong>i7 index</strong> is an 8 base pair barcode specific to a given sample library. <strong>tag</strong> is a 6 base pair tag specific to a cell, and <strong>insert</strong> is the target DNA or RNA fragment from a given sample library.</p>
</div>

Base calling with bc2fastq will produce 3 files per lane: `L001_R1_001.fastq.gz` containing a 6 base pair long cellular barcode, `L001_I1_001.fastq.gz` containing the 8 base pair sample barcode and `L001_R2_001.fastq.gz` containing the reverse complemented 5 prime suffix of the insert region, since it was read in reverse.

>```json
"input": [
    "Lane1_S1_L001_R1_001.fastq.gz",
    "Lane1_S1_L001_I1_001.fastq.gz",
    "Lane1_S1_L001_R2_001.fastq.gz"
],
```
>**declaring input read segments** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling.
{: .example}

Pheniqs could theoretically perform both sample and cellular classification in one run but that can make estimating the priors more tricky. Our downstream pipeline also expected reads for each cell of each sample to be in a separate fastq file. While this may not be ideal for performance and retaining metadata, it is a real world constraint we sometime have to live with to reuse existing workflows we have in place and demonstrates some of  Pheniq's flexibly.

# Phase one: sample classification
In the first phase we will estimate the priors for the sample barcode and split the reads of each sample to a bam file containing the two remaining segments: the cellular barcode and the biological sequence.

>```json
"transform": { "token": [ "0::6", "2::" ] }
```
>**declaring output read segments** Classifying the reads by the sample barcode. We retain the first segment containing the cell tag and the third segment containing the biological sequence.
{: .example}

To classify the reads by the 8 base pair sample barcode we declare the list of possible barcode sequences and a transform that tells Pheniqs where to find the barcode sequence.
>```json
"multiplex": {
  "transform": { "token": [ "1::8" ] }
  "codec": {
    "@GGACTCCT": { "barcode": [ "GGACTCCT" ], "LB": "COL01_PAG069_V2_E2", "SM": "COL01_PAG069" },
    "@ACTCGCTA": { "barcode": [ "ACTCGCTA" ], "LB": "COL20_PAG123_V2_E2", "SM": "COL20_PAG123" },
  }
}
```
>**declaring sample demultiplexing** The transform tells Pheniqs where to locate the sequence that should match the barcode.
{: .example}

To discard reads that failed the internal Illumina sequencer chastity filter we instruct Pheniqs to filter incoming **QC fail** reads
>```json
"filter incoming qc fail": true
```

Since we expect **5%** of reads to be foreign DNA that should not classify to any of the barcodes, we specify the noise prior with the `"noise": 0.05` directive. To mark as **QC fail** read where the posterior probability of correctly decoding the barcode is bellow **0.95** we specify `"confidence threshold": 0.95`.

>```json
{
    "CN": "CGSB AD",
    "PL": "ILLUMINA",
    "PM": "HiSeq",
    "base input url": "~/CBJLFACXX/raw",
    "filter incoming qc fail": true,
    "flowcell id": "CBJLFACXX",
    "flowcell lane number": 1,
    "input": [
        "Lane1_S1_L001_R1_001.fastq.gz",
        "Lane1_S1_L001_I1_001.fastq.gz",
        "Lane1_S1_L001_R2_001.fastq.gz"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "codec": {
            "@COL01_PAG069_V2_E2": {
                "LB": "COL01_PAG069_V2_E2",
                "SM": "COL01_PAG069",
                "barcode": [ "GGACTCCT" ]
            },
            "@COL02_PAG069_V2_E2": {
                "LB": "COL02_PAG069_V2_E2",
                "SM": "COL02_PAG069",
                "barcode": [ "TAAGGCGA" ]
            }
        },
        "confidence threshold": 0.95,
        "noise": 0.05,
        "transform": { "token": [ "1::8" ] }
    },
    "transform": { "token": [ "0::6", "2::" ] }
}
```
>**Phase one** Classifying the sample barcodes using the 3 fastq files produced by bcl2fastq. Here we keep the cellular tag in a separate segment for further processing in the next phase. The actual configuration has more sample barcodes but we only show two here for brevity. The [complete configuration](example/fluidigm/CBJLFACXX_l01_sample.json) is also available.
{: .example}

Before we proceed we validate the configuration.

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json --validate
```

The [validation output](example/fluidigm/CBJLFACXX_l01_sample_validation.txt) is a readable description of all the explicit and implicit parameters after applying defaults. You can check how Pheniqs detects the input format, compression and layout as well as the output you can expect. In The *Output transform* section is a verbal description of how the output read segments will be assembled from the input. similarly *Transform* in the *Mutliplex decoding* section describes how the segment that will be matched against the barcodes is assembled. the You can also see how each of the read groups will be tagged and the prior probability PAMLD will assume for each barcode. You can see that Pheniqs computes the uniform prior.

While not strictly necessary, You may examine the [compile configuration](example/fluidigm/CBJLFACXX_l01_sample_compiled.json), which is the actual configuration Pheniqs will execute with all implicit and default parameters. This is an easy way to see exactly what Pheniqs will be doing and spotting any configuration errors.

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json --compile
```

## Prior estimation

Better estimation of the prior distribution of the samples will improve accuracy. Pheniqs provides a simple python script for adjusting your configuration to include priors estimated from the report emitted by a preliminary Pheniqs run. The `estimate_prior.py` script, distributed with Pheniqs, will either execute Pheniqs with a slightly modified configuration, optimized to refrains from writing the output reads and save time, and emit a modified configuration file with adjusted priors. You may alternatively provide the script a [report from a run you execute yourself]([compile configuration](example/fluidigm/CBJLFACXX_l01_sample_report.json). The priors you specify in your initial configuration can be your best guess for the priors but you can simply leave them out altogether.

>```shell
estimate_prior.py --configuration CBJLFACXX_l01_sample.json \
--split-bam \
--prefix CBJLFACXX_l01
```

or if you already have a report from a preliminary run

>```shell
estimate_prior.py --configuration CBJLFACXX_l01_sample.json \
--report CBJLFACXX_l01_sample_report.json \
--split-bam \
--prefix CBJLFACXX_l01
```

This will produce a [new configuration](example/fluidigm/CBJLFACXX_l01_sample_estimated.json) file with the adjusted priors. For connivence, the script allows you to generate file names in the barcode decelerations so that reads from each barcode are written to a separate file. The method for estimating the prior is [described in the manual](manual.html#prior-estimation).

You can now proceed to demultiplex the sample barcode with the adjusted configuration and produce a separate bam file with reads from each of the samples.

>```shell
pheniqs mux --config CBJLFACXX_l01_sample_estimated.json
```

# Phase two: cellular classification

Our existing pipeline required to have reads from each cell in each sample in a separate fast files. Pheniqs 2.0 does not yet support splitting reads into separate files using cellular barcodes but this is only a minor inconvenience. Since the methodology for decoding the barcodes is essentially the same and since fastq files do not actually contain metadata we can simply use the generic `multiplex` directive to classify the reads by cellular barcodes, in the same fashion we did in the first step. In the first step we produced multiple bam files, one for each sample. In this second phase we will estimate priors and decode cellular barcodes on each of those bam files independently so the following should be applied to each of the bam files produced in the first phase.

First we notice that this time we have 2 input segments. The first is the 6 base pair cellular barcode and the second is our DNA or RNA fragment.

>```json
"transform": { "token": [ 1::" ] }
```
>**declaring output read segments** Classifying the reads by the cellular barcode. We retain the second segment containing the biological sequence.
{: .example}

To classify the reads by the 6 base pair sample barcode we declare the list of possible barcode sequences and a transform that tells Pheniqs where to find the barcode sequence.

>```json
"multiplex": {
  "transform": { "token": [ "0::6" ] }
  "codec": {
    "@CACGTA": { "barcode": [ "CACGTA" ]},
    "@CTCACA": { "barcode": [ "CTCACA" ]}
  }
}
```
>**declaring cellular demultiplexing** The transform tells Pheniqs where to locate the sequence that should match the barcode.
{: .example}

Since the cellular barcodes are identical for all first phase bam files we intentionally refrain from declaring the input directive in the configuration file and will provide it on the command line so that we can reuse the same configuration file all, sample specific, bam files.

>```json
{
    "CN": "CGSB AD",
    "PL": "ILLUMINA",
    "PM": "HiSeq",
    "base input url": "~/CBJLFACXX/sample",
    "base output url": "~/CBJLFACXX/cellular",
    "transform": { "token": [ "1::" ] },
    "flowcell id": "CBJLFACXX",
    "filter incoming qc fail": true,
    "filter outgoing qc fail": true,
    "multiplex": {
        "algorithm": "pamld",
        "confidence threshold": 0.95,
        "noise": 0.05,
        "transform": { "token": [ "0::6" ] },
        "codec": {
            "@CACGTA": { "barcode": [ "CACGTA" ]},
            "@CTCACA": { "barcode": [ "CTCACA" ]}
        }
    }
}
```
>**Phase two** Classifying the cellular barcodes. The actual configuration has more sample barcodes but we only show two here for brevity. The [complete configuration](example/fluidigm/CBJLFACXX_l01_cellular.json) is also available.
{: .example}

An advantage of decoding cellular barcodes on each sample separately is that we can estimate the sample conditional priors.

>```shell
estimate_prior.py --configuration CBJLFACXX_l01_cellular.json \
--input CBJLFACXX_l01_COL01_PAG069_V2_E2.bam \
--prefix CBJLFACXX_l01_COL01_PAG069_V2_E2 \
--split-fastq \
--sense-input
```
>**estimating cellular priors** We execute this for each of the bam files produced in the first phase. the `--input` command line parameter for `estimate_prior.py` will be forwarded to Pheniqs when estimating the priors and will also be present in the generated configuration with the adjusted priors. The `--split-fastq` parameter will declare an `output` directive in each of the barcodes so that the configuration produce the desired, split fastq, output.
{: .example}

If you already have a [report](example/fluidigm/CBJLFACXX_l01_COL01_PAG069_V2_E2_report.json) from a preliminary run.

>```shell
estimate_prior.py --configuration CBJLFACXX_l01_cellular.json \
--report CBJLFACXX_l01_COL01_PAG069_V2_E2_report.json \
--input CBJLFACXX_l01_COL01_PAG069_V2_E2.bam CBJLFACXX_l01_COL01_PAG069_V2_E2.bam \
--prefix CBJLFACXX_l01_COL01_PAG069_V2_E2 \
--split-fastq
```
>**estimating cellular priors with existing report** Since the `--input` command line parameter is forwarded to Pheniqs we specify the input twice to tell Pheniqs we expect two segments to be pulled from the file. Using the automatic input sensing can be simple but requires access to the actual input file.
{: .example}

This will produce an [adjusted configuration](example/fluidigm/CBJLFACXX_l01_COL01_PAG069_V2_E2_adjusted.json) with the estimated priors, specific to this particular sample.


>```json
{
    "CN": "CGSB AD",
    "PL": "ILLUMINA",
    "PM": "HiSeq",
    "base input url": "~/CBJLFACXX/sample",
    "base output url": "~/CBJLFACXX/cellular",
    "filter incoming qc fail": true,
    "filter outgoing qc fail": true,
    "flowcell id": "CBJLFACXX",
    "input": [
        "CBJLFACXX_l01_COL01_PAG069_V2_E2.bam?format=bam",
        "CBJLFACXX_l01_COL01_PAG069_V2_E2.bam?format=bam"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "codec": {
            "@ACACTG": {
                "barcode": [
                    "ACACTG"
                ],
                "concentration": 0.04284168923365803,
                "output": [
                    "CBJLFACXX_l01_COL01_PAG069_V2_E2_ACACTG_s01.fastq.gz"
                ]
            },
            "@ACATGC": {
                "barcode": [
                    "ACATGC"
                ],
                "concentration": 0.031949147727610774,
                "output": [
                    "CBJLFACXX_l01_COL01_PAG069_V2_E2_ACATGC_s01.fastq.gz"
                ]
            }
        },
        "confidence threshold": 0.95,
        "noise": 0.0990339866699282,
        "transform": {
            "token": [
                "0::6"
            ]
        }
    },
    "transform": {
        "token": [
            "1::"
        ]
    }
}
```
>**Classifying cellular barcodes with estimated priors** The actual configuration has more cellular barcodes but we only show two here for brevity. The [complete configuration](example/fluidigm/CBJLFACXX_l01_COL01_PAG069_V2_E2_adjusted.json) is also available.
{: .example}

Executing Pheniqs with this configuration will generate the desired individual fastq files, each contaning only reads for one cellular barcode in one sample.
