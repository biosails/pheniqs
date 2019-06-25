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

## Basecalling

The `pheniqs-illumina-recipe.py` script can help you execute bcl2fastq to perform basecalling without barcode decoding.

```
pheniqs-illumina-recipe.py basecall --fastq-compression-level 3 181014_A00534_0024_AH7LT2DSXX
```

will write `H7LT2DSXX_basecall.sh` and `basecall_samplesheet.csv` to the current folder. You can also provide values for some relevant *bcl2fastq* parameters.

>```shell
    bcl2fastq \
    --runfolder-dir \
    illumina/181014_A00534_0024_AH7LT2DSXX \
    --sample-sheet basecall_samplesheet.csv \
    --create-fastq-for-index-reads \
    --adapter-stringency 0 \
    --minimum-trimmed-read-length 0 \
    --mask-short-adapter-reads 0 \
    --fastq-compression-level 3
```
>**example bcl2fastq basecall shell script** We must provide bcl2fastq an alternative sample sheet or it will default to using the one in the run folder. A simple sample sheet that does not perform any barcode decoding [is generated]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/basecall_samplesheet.csv) by `pheniqs-illumina-recipe.py basecall`. You may also choose the gzip compression level. For temporary files it is better to choose low values since IO will be faster and file size only marginally larger.
{: .example}

## Core configuration

[H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) is a core configuration file that summarizes metadata extracted from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml) and [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv). It is imported by most other configuration files in this tutorial. To generate one use the `core` subcommand of `pheniqs-illumina-api.py`

>```shell
pheniqs-illumina-api.py core illumina/181014_A00534_0024_AH7LT2DSXX --no-input-npf
```
>**Generating a core configuration** `--no-input-npf` will add a global `filter incoming qc fail` instruction to discard reads that failed the internal Illumina sequencer chastity filter from the incoming feed.
{: .example}

*bc2fastq* will produced 4 files per lane:

>```json
"input": [
  "H7LT2DSXX_S1_L001_R1_001.fastq.gz",
  "H7LT2DSXX_S1_L001_I1_001.fastq.gz",
  "H7LT2DSXX_S1_L001_I2_001.fastq.gz",
  "H7LT2DSXX_S1_L001_R2_001.fastq.gz"
]
```
>**input read segments** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling. `R1` containing the 3 prime prefix of the insert region, `I1` containing the i7 index, `I2` containing the i5 index and `R2` containing the reverse complemented 5 prime suffix of the insert region, since it was read in reverse.
{: .example}

To classify the reads by the i5 and i7 indices we need to declare the list of possible barcode sequences and a tokens that tell Pheniqs where to find the barcode sequence. The [core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) configuration file conveniently declares an array of decoders that where recovered from the [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv) file as well as transformation rules from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml). If you import `H7LT2DSXX_core.json` You can reuse those in your configuration file and expand them so you don't need to constantly be editing big configuration files.

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

In `H7LT2DSXX_core.json` you will find a `decoder` directive that defines a dictionary of abstract decoders.

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
>**H7LT2DSXX_l01_multiplex decoder** lists the possible barcode combinations and the library names associated with them in the [LB](glossary.html#lb_auxiliary_tag) tag extracted from
[SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv).
{: .example}

You can use those decoders as a starting point in the `multiplex`, `cellular`, and `molecular` directive with `base`. Any directive you specify in your instantiation will override values provided by the base.

## Decoding without priors

The `multiplex` sub command will generate a basic sample demultiplexing configuration file for each lane.

```
pheniqs-illumina-api.py multiplex illumina/181014_A00534_0024_AH7LT2DSXX
```

To decode sample barcodes without prior estimation we declare a `multiplex` directive that expands the `H7LT2DSXX_l01_multiplex` decoder that we have seen defined in `H7LT2DSXX_core.json`. The report is written to a file instead of standard error with `report url`.  You should still specify your best guess for the `noise` prior. For example [H7LT2DSXX_l01_sample.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_sample.json) is a uniform configuration for the first lane.

>```json
{
    "import": [
        "H7LT2DSXX_core.json"
    ],
    "input": [
        "H7LT2DSXX_S1_L001_R1_001.fastq.gz",
        "H7LT2DSXX_S1_L001_I1_001.fastq.gz",
        "H7LT2DSXX_S1_L001_I2_001.fastq.gz",
        "H7LT2DSXX_S1_L001_R2_001.fastq.gz"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "base": "H7LT2DSXX_l01_multiplex",
        "confidence threshold": 0.95,
        "noise": 0.05
    },
    "output": [
        "H7LT2DSXX_l01.bam"
    ],
    "report url": "H7LT2DSXX_l01_sample_report.json",
    "transform": {
        "token": [
            "0::",
            "3::"
        ]
    }
}```
>**Decoding with a uniform prior** output is written to a bam file.
{: .example}

## Prior estimation

The `estimate` sub command will generate a sample prior estimation optimized configuration file for each lane.

```
pheniqs-illumina-api.py estimate illumina/181014_A00534_0024_AH7LT2DSXX
```

To estimate priors for sample barcodes we need to collect statistics about the sequences identified by the tokens. We don't actually need to read the 2 segments containing the biological sequences so we declare input only for the two segments containing the indices. We then declare a `multiplex` directive that expands the `H7LT2DSXX_l01_multiplex` decoder that we have seen defined in `H7LT2DSXX_core.json`. The correct tokenization for the modified input is the first 8 bases of segment 0 and 1. `output` is redirected to `/dev/null` to tell Pheniqs it should not bother with the output. The report is written to a file instead of standard error with `report url`. For example [H7LT2DSXX_l01_estimate.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_estimate.json) is a uniform configuration for the first lane.

>```json
{
    "import": [
        "H7LT2DSXX_core.json"
    ],
    "input": [
        "H7LT2DSXX_S1_L001_I1_001.fastq.gz",
        "H7LT2DSXX_S1_L001_I2_001.fastq.gz"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "base": "H7LT2DSXX_l01_multiplex",
        "confidence threshold": 0.95,
        "noise": 0.05,
        "transform": {
            "token": [
                "0::",
                "1::"
            ]
        }
    },
    "output": [
        "/dev/null"
    ],
    "report url": "H7LT2DSXX_l01_estimate_report.json",
    "transform": {
        "token": [
            "0::",
            "1::"
        ]
    }
}
```
>**Prior estimation configuration** refraining from reading the biological sequences and producing no output significantly speeds things up.
{: .example}

Executing this configuration will yield the report
[H7LT2DSXX_l01_estimate_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/sample/prior/H7LT2DSXX_l01_estimate_report.json). The [standard Pheniqs report](manual.html#quality-control-and-statistics) contains decoding statistics that are used to estimate the priors.

>```shell
pheniqs mux --config sample/prior/l01_sample.json
```

The report can now be used to generate an adjusted configuration from the uniform one

>```shell
estimate_prior.py \
--report H7LT2DSXX_l01_estimate_report.json \
--configuration sH7LT2DSXX_l01_sample.json \
> H7LT2DSXX_l01_adjusted.json
```
>**decoding a lane with an estimated prior** [H7LT2DSXX_l01_adjusted.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_adjusted.json) is similar to [H7LT2DSXX_l01_sample.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_sample.json) with the addition of the estimated priors.
{: .example}

## Decoding with the estimated prior

>```shell
pheniqs mux --config H7LT2DSXX_l01_adjusted.json
```
>**decoding a lane with an estimated prior** the report will be written to [H7LT2DSXX_l01_adjusted_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_adjusted_report.json) because it is relative to the base output.
{: .example}
