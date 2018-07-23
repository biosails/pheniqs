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
        <li><a class="active"   href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/2.0/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Quick Tutorial
{:.page-title}

* placeholder
{:toc }

The Pheniqs command line interface accepts a [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration file containing a number of separate sections specifying directives for input and output layout, parsing read segments and run parameters. In the [workflow page](workflow.html) you will find some annotated examples of complete configuration files. If you are new to JSON, a [validator](cli.html#json-validation) can be instrumental for troubleshooting syntax errors. Some parameters are also exposed as command line arguments that override their corresponding configuration file values. A brief description of the command line parameters Pheniqs accepts is always available with the `-h/--help` flags. If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to [install the bundled command line completion](cli.html#zsh-completion) script for a more interactive command line experience.

# File Format
Pheniqs can arbitrarily manipulate reads from either [SAM, BAM and CRAM](glossary.html#htslib) or [FASTQ](glossary.html#fastq) with segments either [interleaved](glossary.html#interleaved_file_layout) into a single file or [split](glossary.html#split_file_layout) over many. Read manipulation is achieved by means of [tokenization](#tokenization) and [construction](#construction). In the tokenization step Pheniqs consults the token patterns you declared to extract tokens from an [input segment](glossary.html#input_segment). In the construction step [transform patterns](manual.html#transform-pattern) reference the token patterns to construct new segments. The optional construction directive is only necessary when composing output segments from multiple, non continuous, tokens and if omitted each token is assumed to declare a single output segment.

## The input directive
The simplest configuration can include just an `input` directive. In this example we consider three files that contain synchronized segments from an Illumina MiSEQ instrument. Pheniqs will assemble an input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input feed.

>```json
{
    "input": [
        "BDGGG_s01.fastq.gz",
        "BDGGG_s02.fastq.gz",
        "BDGGG_s03.fastq.gz"
    ]
}
```
>Declaring an input read that is [split](glossary.html#split_file_layout) over three gzip compressed FASTQ files.
{: .example}

This is already a valid Pheniqs configuration! Since no manipulation instructions are specified reads are simply interleaved to the output. Since output is not explicitly declared it defaults to the SAM format and written to standard output.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     FI:i:1  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   13      *       0       0       *       *       0       0       GGACTCCT        B@CCCFC<        FI:i:2  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     FI:i:3  TC:i:3
```
>Output from [interleaving](glossary.html#interleaved_file_layout) the three read segments into a single unmanipulated SAM formatted stream written to standard output.
{: .example}

This simple configuration can be very useful for packing raw split read segments into a single CRAM file. CRAM files are the latest indexed and compressed binary encoding of SAM and are much faster to interact with. Packing your reads into CRAM also makes archiving easier and will often use up less storage space. To write the interleaved output to a compressed CRAM file simply add an output directive.

>```json
{
    "input": [
        "BDGGG_s01.fastq.gz",
        "BDGGG_s02.fastq.gz",
        "BDGGG_s03.fastq.gz"
    ],
    "output": [ "BDGGG_raw.cram" ]
}
```
>Interleaving three unmanipulated read segments into a single efficient CRAM file. CRAM files are often much faster to read and write, especially in highly parallelized environments, and also support a rich metadata vocabulary.
{: .example}

## Tokenization
The `template` directive can be used to manipulate the structure of the output read. If omitted all segments of the input are written verbatim to the output, which is what we did when we interleaved our 3 FASTQ files into a single CRAM file. Since the second segment contains only a technical sequence and we don't wish to write it to the output we add instructions to assemble an output read from only the first and third segments of the input.

>```json
{
    "input": [
        "BDGGG_s01.fastq.gz",
        "BDGGG_s02.fastq.gz",
        "BDGGG_s03.fastq.gz"
    ],
    "template": { "token": [ "0::", "2::" ] }
}
```
>Adding a template directive composing an output read from only the untouched first and third input segments. Input segments in pheniqs are indexed and referenced using a [zero based coordinate system](glossary.html#zero_based_coordinate) so the first segment is 0.

The [token patterns](manual.html#tokenization) declared in the token array of the template directive are made of 3 colon separated integers. The first is the [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment). The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. The two colons are always mandatory. Pheniqs token pattern can address segments from either the 5' (left) or 3' (right) end and mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
```
>Output from interleaving the three read segments to a SAM formatted stream. Notice that only the first and third segments were written.
{: .example}

## Read Group classification
The reads in our files are DNA sequenced from 5 different libraries that were tagged with an 8bp technical sequence and pooled together for sequencing. We want to classify the reads into read groups by examining the first 8 nucleotides of the second segment. Decoding the barcode would be a trivial matter of comparing two character sequences if we were absolutely confident no errors were made during sequencing. Sadly this is not the case and each nucleotide reported by a sequencing instrument is accompanied by an estimate of the probability the base was incorrectly called. In the `multiplex` directive we declare a single decoder used to classify the reads into read groups.

>```json
{
    "multiplex": {
        "template": { "token": [ "1::8" ] },
        "codec": {
            "@AGGCAGAA": { "barcode": [ "AGGCAGAA" ] },
            "@CGTACTAG": { "barcode": [ "CGTACTAG" ] },
            "@GGACTCCT": { "barcode": [ "GGACTCCT" ] },
            "@TAAGGCGA": { "barcode": [ "TAAGGCGA" ] },
            "@TCCTGAGC": { "barcode": [ "TCCTGAGC" ] }
        }
        "algorithm": "pamld",
        "noise": 0.02,
        "confidence threshold": 0.95,
    }
}
```
>A `multiplex` decoder directive using the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The embedded `template` directive is used to extract the barcode from the raw read segments specifying while the `codec` directive names the possible barcode sequences we expect to find.

In this example we declare a `multiplex` directive that uses the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) algorithm. This algorithm will reject any classification with a decoding confidence lower than the threshold specified in `confidence threshold`. The `noise` parameter is the prior probability that a sequence was not derived from any of the provided barcodes and is often set to the amount of [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution for the Illumina platform. In the `codec` dictionary we provide a discrete set of possible decoding results. All `barcode` arrays must match the layout declared in the embedded `template` directive, in our case one 8bp segment. Degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) bases are currently not allowed in a multiplex barcode sequence. *The keys of the `codec` dictionary can be any unique string and in this example we specified them to be the unique barcode nucleotide sequence prefixed with an @ character to remind you this is simply a unique identifier.*

>```json
{
    "input": [
        "BDGGG_s01.fastq.gz",
        "BDGGG_s02.fastq.gz",
        "BDGGG_s03.fastq.gz"
    ],
    "template": { "token": [ "0::", "2::" ] },
    "multiplex": {
        "algorithm": "pamld",
        "noise": 0.02,
        "confidence threshold": 0.95,
        "template": { "token": [ "1::8" ] },
        "codec": {
            "@AGGCAGAA": { "barcode": [ "AGGCAGAA" ], "LB": "trinidad 5" },
            "@CGTACTAG": { "barcode": [ "CGTACTAG" ], "LB": "trinidad 4" },
            "@GGACTCCT": { "barcode": [ "GGACTCCT" ], "LB": "trinidad 9" },
            "@TAAGGCGA": { "barcode": [ "TAAGGCGA" ], "LB": "trinidad 1" },
            "@TCCTGAGC": { "barcode": [ "TCCTGAGC" ], "LB": "trinidad 8" }
        }
    },
    "CN": "CGSB",
    "DT": "2018-02-25T07:00:00+00:00",
    "PI": "300",
    "PL": "ILLUMINA",
    "PM": "miseq",
    "SM": "trinidad",
    "flowcell id": "BDGGG",
    "flowcell lane number": 1
}
```
> Demultiplexing with one 8bp barcode on the second segment to an [interleaved](glossary.html#interleaved_file_layout) CRAM file. The [CN](glossary.html#cn_auxiliary_tag), [DT](glossary.html#dt_auxiliary_tag), [PI](glossary.html#pi_auxiliary_tag), [PL](glossary.html#pl_auxiliary_tag), [PM](glossary.html#pm_auxiliary_tag) and [SM](glossary.html#sm_auxiliary_tag) tags are declared globally and will be added to all read groups. The [LB](glossary.html#lb_auxiliary_tag) is declared individually for each barcode.
{: .example}

The SAM [RG](glossary.html#rg_auxiliary_tag) header tag can contain additional metadata. Capitalized 2 letter directives in pheniqs configuration files often refer to their corresponding SAM tag and can be either specified with the individual barcodes or globally for inclusion in all of them.


>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:BDGGG:1:AGGCAGAA     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 5   PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:AGGCAGAA     SM:trinidad
@RG     ID:BDGGG:1:CGTACTAG     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 4   PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:CGTACTAG     SM:trinidad
@RG     ID:BDGGG:1:GGACTCCT     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 9   PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:GGACTCCT     SM:trinidad
@RG     ID:BDGGG:1:TAAGGCGA     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 1   PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:TAAGGCGA     SM:trinidad
@RG     ID:BDGGG:1:TCCTGAGC     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 8   PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:TCCTGAGC     SM:trinidad
@RG     ID:BDGGG:1:undetermined CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:1:undetermined SM:trinidad
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     RG:Z:BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     RG:Z:BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
M02455:162:000000000-BDGGG:1:1101:10000:12232   77      *       0       0       *       *       0       0       GTATAGGGGTCACATATAGTTGGTGTGCTTTGTGAACTGCGATCTTGACGG     CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     RG:Z:BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:CCCCCGGG   XB:f:1.56086e-06
M02455:162:000000000-BDGGG:1:1101:10000:12232   141     *       0       0       *       *       0       0       GTCCTATCCTACTCGGCTTCTCCCCATTTTTCAGACATTTTCCTATCAGTC     CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     RG:Z:BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:CCCCCGGG   XB:f:1.56086e-06
```
>Output from demultiplexing using the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) using the above configuration. Notice how tags declared globally were added to all read groups. The [RG](glossary.html#rg_auxiliary_tag) and [PU](glossary.html#pu_auxiliary_tag) read group identifiers default to the convention set by [GATK](https://software.broadinstitute.org/gatk/guide/article?id=6472) if you provide the flowcell related directives.
{: .example}

# Decoding
Pheniqs implements a novel [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that directly estimates the decoding likelihood from the base calling error probabilities provided by the sequencing platform. When output is SAM encoded Pheniqs will write the multiplex barcode decoding error probability to the [XB](glossary.html#dq_auxiliary_tag) auxiliary tag. Pheniqs also implements a traditional [minimum distance decoder](glossary.html#minimum_distance_decoding) that consults the edit distance between the expected and observed sequence. To select a decoding algorithm set `algorithm` property to either `pamld` or `mdd`.

When using phred-adjusted maximum likelihood decoding you may want to set the `confidence threshold` and `noise` properties. The `confidence threshold` property is the threshold on the decoding probability for the decoder to declare a successful classification. If the decoder fails to classify the read it is considered undetermined. The value of `confidence threshold` defaults to **0.99** but in practice depends on the application. Since the multiplex barcode decoding error probability is written to the [XB](glossary.html#dq_auxiliary_tag) auxiliary tag you can set `confidence threshold` to a relatively lax value and leave the decision to exclude the read for downstream analysis.

The `noise` property can be used to provide a prior probability of observing an undetermined read. For the Illumina platform the most common value will be the concentration of the [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution, which is usually 1% and corresponds to a value of **0.01**. The `noise` property defaults to **0.01**.
