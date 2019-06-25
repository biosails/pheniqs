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

This tutorial will walk you through prior estimation and demultiplexing of an Illumia NovaSeq 6000 sequencing run with [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The run has paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol. Each read has 4 segments: two 151 base pairs long biological sequences from the DNA or RNA fragment and two 8 base pairs long technical sequences containing the i5 and i7 indices.

![paird end sequencing](/pheniqs/assets/img/paired_end_sequencing.png)

In this example we will be using `pheniqs-illumina-api.py` to generate configuration files from metadata present in the [Illumina run folder]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX). We will also be demonstrating using the `import` directive for reusing configuration instructions. For each read, the sample barcode and its quality scores will be written to the [BC tag](glossary.html#bc_auxiliary_tag) and [QT tag](glossary.html#qt_auxiliary_tag) respectively, and the decoding error probability to the [XB tag](glossary.html#xb_auxiliary_tag). **All shell commands bellow are executed in the [H7LT2DSXX example directory]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX), where you can also find pre generated files for this example**.

## Basecalling

Generate a bcl2fastq shell command for base calling without decoding the sample barcodes with `pheniqs-illumina-api.py`.

>```shell
    pheniqs-illumina-api.py basecall \
    --output-dir . \
    --fastq-compression-level 3 \
    181014_A00534_0024_AH7LT2DSXX
```
>**generating a bcl2fastq command** we want output fastq files to be written to this directory and we want lower gzip compression since we will not be archiving those files.
{: .example}

This will produce [H7LT2DSXX_basecall.sh]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_basecall.sh) and [H7LT2DSXX_basecall_sample_sheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/HH7LT2DSXX_basecall_sample_sheet.csv). `pheniqs-illumina-api.py basecall` allows you to forward several relevant *bcl2fastq* parameters to be included in the script, but you can also manually adjust it.

>```shell
    bcl2fastq \
    --runfolder-dir 181014_A00534_0024_AH7LT2DSXX \
    --sample-sheet H7LT2DSXX_basecall_sample_sheet.csv \
    --create-fastq-for-index-reads \
    --adapter-stringency 0 \
    --minimum-trimmed-read-length 0 \
    --mask-short-adapter-reads 0 \
    --output-dir . \
    --fastq-compression-level 3

```
>**bcl2fastq basecall shell script** We must provide bcl2fastq an alternative sample sheet or it will default to using the one in the run folder. A simple sample sheet that does not perform any barcode decoding [is generated]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/HH7LT2DSXX_basecall_sample_sheet.csv) by `pheniqs-illumina-api.py basecall`.
{: .example}

Base calling with bcl2fastq took about 6 hours on a *dual socket Intel Xeon E5-2620* with a total of 12 cores and produced 4 files per lane.

>```shell
H7LT2DSXX_S1_L001_I1_001.fastq.gz
H7LT2DSXX_S1_L001_I2_001.fastq.gz
H7LT2DSXX_S1_L001_R1_001.fastq.gz
H7LT2DSXX_S1_L001_R2_001.fastq.gz
H7LT2DSXX_S1_L002_I1_001.fastq.gz
H7LT2DSXX_S1_L002_I2_001.fastq.gz
H7LT2DSXX_S1_L002_R1_001.fastq.gz
H7LT2DSXX_S1_L002_R2_001.fastq.gz
H7LT2DSXX_S1_L003_I1_001.fastq.gz
H7LT2DSXX_S1_L003_I2_001.fastq.gz
H7LT2DSXX_S1_L003_R1_001.fastq.gz
H7LT2DSXX_S1_L003_R2_001.fastq.gz
```
>**bcl2fastq output** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling. `R1` containing the 3 prime prefix of the insert region, `I1` containing the i7 index, `I2` containing the i5 index and `R2` containing the reverse complemented 5 prime suffix of the insert region, since it was read in reverse.
{: .example}

## Core configuration

[H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) summarizes metadata extracted from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml) and [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv). It is imported by most other configuration files in this tutorial. To generate one use the `core` subcommand of `pheniqs-illumina-api.py`

>```shell
pheniqs-illumina-api.py core \
--base-input-url . \
--base-output-url . \
--no-input-npf \
181014_A00534_0024_AH7LT2DSXX
```
>**Generating a core configuration** `--no-input-npf` will add a global `filter incoming qc fail` instruction to discard reads that failed the internal Illumina sequencer chastity filter from the incoming feed.
{: .example}

[H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) conveniently declares an array of decoders that where recovered from the [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/illumina/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv) file as well as transformation rules from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml). When you import [H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) you can reuse those in your configuration file and expand them so you don't need to constantly be editing big configuration files. Anything specified in a configuration file will override an imported directive.

>```json
{
    "PL": "ILLUMINA",
    "PM": "A00534",
    "flowcell id": "H7LT2DSXX"
}
```
>**core configuration** already declares [PL](glossary.html#pl_auxiliary_tag),
[PM](glossary.html#pm_auxiliary_tag), and `flowcell id`. Those were extracted from [RunInfo.xml]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX/RunInfo.xml) but you can modify them.
{: .example}

## decoder configuration element

[H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json) declares a `decoder` directive that defines a dictionary of abstract decoders.

>```json
{
    "decoder": {
        "H7LT2DSXX_l01_multiplex": {
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
>**H7LT2DSXX_l01_multiplex decoder** lists the possible barcode combinations and the library names associated with them in the [LB](glossary.html#lb_auxiliary_tag) tag extracted from [SampleSheet.csv]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv).
{: .example}

You can the `base` property on a decoder to use those decoders as a starting point in the `multiplex`, `cellular`, and `molecular`. Any directive you specify in your instantiation will override values provided by the referenced base.

## Decoding without priors

The `multiplex` sub command will generate a basic sample demultiplexing configuration file for each lane.

>```shell
pheniqs-illumina-api.py multiplex \
--confidence 0.95 \
--noise 0.05 \
181014_A00534_0024_AH7LT2DSXX
```
>**pheniqs-illumina-api multiplex** will produce a basic sample demultiplexing configuration that we will adjust later with the prior estimates. Some pheniqs parameters can be forwarded like `--noise` and `--confidence` shown in this example.
{: .example}

[H7LT2DSXX_l01_sample.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_sample.json)
is a configuration for demuliplexing the first lane. It declares a `multiplex` directive that expands the `H7LT2DSXX_l01_multiplex` decoder we have seen defined in [H7LT2DSXX_core.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_core.json).

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
    "report url": "H7LT2DSXX_l01_sample_report.json"
}
```
>**Decoding with a uniform prior** output is written to a bam file.
{: .example}

## Prior estimation

The `estimate` sub command will generate a sample prior estimation optimized configuration file for each lane.

>```shell
pheniqs-illumina-api.py estimate \
--confidence 0.95 \
--noise 0.05 \
181014_A00534_0024_AH7LT2DSXX
```
>**pheniqs-illumina-api estimate** will produce a sample demultiplexing configuration optimized for prior estimation.
{: .example}

To estimate priors for sample barcodes we need to collect some statistics about the sample barcodes. Since don't need to read the 2 segments containing the biological sequences input is declared only for the two segments containing the indices.

[H7LT2DSXX_l01_estimate.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_estimate.json) declares a `multiplex` directive that expands `H7LT2DSXX_l01_multiplex` and adjusts the tokenization for the modified input. `output` is redirected to `/dev/null` to tell Pheniqs it should not bother with the output.

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
>**Prior estimation configuration** refrains from reading the biological sequences and produces no output which significantly speeds things up..
{: .example}

Executing this configuration will yield the report
[H7LT2DSXX_l01_estimate_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_estimate_report.json). Like every [Pheniqs report](manual.html#quality-control-and-statistics), it contains decoding statistics that we can use to estimate the priors.

>```shell
pheniqs mux --config H7LT2DSXX_l01_estimate.json
```
>**executing pheniqs with a prior estimation configuration** will produce [H7LT2DSXX_l01_estimate_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_estimate_report.json).
{: .example}

[H7LT2DSXX_l01_estimate_report.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_estimate_report.json) and [H7LT2DSXX_l01_sample.json]({{ site.github.repository_url }}/blob/master/example/H7LT2DSXX/H7LT2DSXX_l01_sample.json) can now be used to generate an adjusted configuration.

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
