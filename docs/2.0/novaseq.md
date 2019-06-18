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
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Build</a></li>
        <li><a class="active"   href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/best_practices.html">Best Practice</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# NovaSeq 6000
{:.page-title}

This tutorial will walk you through demultiplexing an Illumia NovaSeq 600 sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The read has 4 segments: 2x151bp biological from the DNA or RNA fragment and 2x8bp technical containing the i5 and i7 indices. If the results are written to SAM, BAM or CRAM the multiplex barcode and its quality scores are written to the [BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag) respectively, and the decoding error probabilities is written to the [XB](glossary.html#xb_auxiliary_tag).

![paird end sequencing](/pheniqs/assets/img/paired_end_sequencing.png)

## Core configuration

[core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/core.json) is a core configuration file that summarizes metadata extracted from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml) and [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv). It is imported by most other configuration files in this tutorial. To generate one use the `core` subcommand of `illumina2pheniqs.py`

```
illumina2pheniqs.py core illumina/181014_A00534_0024_AH7LT2DSXX > core.json
```

## Input Read Layout

Base calling with bc2fastq produced 4 files per lane:

>```json
"input": [
  "H7LT2DSXX_l01_S1_L001_R1_001.fastq.gz",
  "H7LT2DSXX_l01_S1_L001_I1_001.fastq.gz",
  "H7LT2DSXX_l01_S1_L001_I2_001.fastq.gz",
  "H7LT2DSXX_l01_S1_L001_R2_001.fastq.gz"
]
```
>**input read segments** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling. `R1` containing the 3 prime prefix of the insert region, `I1` containing the i7 index, `I2` containing the i5 index and `R2` containing the reverse complemented 5 prime suffix of the insert region, since it was read in reverse.
{: .example}

To classify the reads by the i5 and i7 indices we need to declare the list of possible barcode sequences and a tokens that tell Pheniqs where to find the barcode sequence. The [core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/core.json) configuration file conveniently declares an array of decoders that where recovered from the [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv) file as well as transformation rules from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml). If you import `core.json` You can reuse those in your configuration file and expand them so you don't need to constantly be editing big configuration files.

>```json
{
    "PL": "ILLUMINA",
    "PM": "A00534",
    "flowcell id": "H7LT2DSXX"
}
```
>**core configuration** already declares [PL](glossary.html#pl_auxiliary_tag),
[PM](glossary.html#pm_auxiliary_tag),
and the `flowcell id`. Those were extracted from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml) but you can modify them.
{: .example}

## decoder configuration element

In `core.json` you will find a `decoder` directive that defines a dictionary of abstract decoders.

>```json
{
    "decoder": {
        "H7LT2DSXX_lane_1_multiplex": {
            "transform": { "token": [ "1::8", "2::8" ] },
            "codec": {
                "@A10_PDAC81": {
                    "LB": "A10_PDAC81",
                    "barcode": [
                        "CGAGGCTG",
                        "GTAAGGAG"
                    ]
                },
                "@A11_PDAC490": {
                    "LB": "A11_PDAC490",
                    "barcode": [
                        "AAGAGGCA",
                        "ACTGCATA"
                    ]
                }
            }
        }
    }
}
```
>**H7LT2DSXX_lane_1_multiplex decoder** lists the possible barcode combinations and the library names associated with them in the [LB](glossary.html#lb_auxiliary_tag) tag extracted from
[SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv).
{: .example}

You can use those decoders as a starting point in the `multiplex`, `cellular`, and `molecular` directive with `base`. Any directive you specify in your instantiation will override values provided by the base.

## Prior estimation

To estimate priors for sample barcodes we need to collect statistics about the sequences identified by the tokens. We don't actually need to read the 2 segments containing the biological sequences so we declare input only for the two segments containing the indices. We then declare a `multiplex` directive that expands the `H7LT2DSXX_lane_1_multiplex` decoder that we have seen defined in `core.json`. The correct tokenization for the two inputs is the first 8 bases of segment 0 and 1. `output` is redirected to `/dev/null` to tell Pheniqs it should not bother with the output. The report is written to a file instead of standard error with `report url`. To discard reads that failed the internal Illumina sequencer chastity filter from the incoming feed we specify `filter incoming qc fail`.

>```json
{
    "filter incoming qc fail": true,
    "report url": "l01_sample_report.json",
    "import": [
        "../../core.json"
    ],
    "input": [
        "H7LT2DSXX_l01_S1_L001_I1_001.fastq.gz",
        "H7LT2DSXX_l01_S1_L001_I2_001.fastq.gz"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "base": "H7LT2DSXX_lane_1_multiplex",
        "confidence threshold": 0.95,
        "noise": 0.05,
        "transform": {
            "token": [
                "0::8",
                "1::8"
            ]
        }
    },
    "output": [
        "/dev/null"
    ]
}
```
>**Prior estimation configuration** refraining from reading the biological sequences and producing no output significantly speeds things up.
{: .example}

Executing this configuration will yield the report
[l01_sample_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/sample/prior/l01_sample_report.json). This is the standard report Pheniqs produces and it contains decoding statistics that are used to estimate the priors.

>```shell
pheniqs mux --config sample/prior/l01_sample.json \
--base-input ~/H7LT2DSXX \
--base-output sample/prior
```

Producing an adjusted configuration from the report

>```shell
estimate_prior.py \
--report sample/prior/l01_sample_report.json \
--configuration sample/uniform/l01_sample.json \
> sample/adjusted/l01_sample.json
```

To emit the two ends of the insert region as two segments of the output read we declare the global transform directives
>```json
"transform": { "token": [ "0::", "3::" ] }
```
>**declaring output read segments** Only the segments coming from the first and the forth file are biological sequences and should be included in the output.
{: .example}
