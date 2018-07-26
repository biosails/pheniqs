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

# Supported File Format
Pheniqs can arbitrarily manipulate reads from either [SAM, BAM and CRAM](glossary.html#htslib) or [FASTQ](glossary.html#fastq) with segments either [interleaved](glossary.html#interleaved_file_layout) into a single file or [split](glossary.html#split_file_layout) over many. Read manipulation is achieved by means of [tokenization](#tokenization) and [construction](#construction). In the tokenization step Pheniqs consults the token patterns you declared to extract tokens from an [input segment](glossary.html#input_segment). In the construction step [transform patterns](manual.html#transform-pattern) reference the token patterns to construct new segments. The optional construction directive is only necessary when composing output segments from multiple, non continuous, tokens and if omitted each token is assumed to declare a single output segment.

# Declaring Input
A very simple configuration can include nothing more than an `input` directive. In this example we consider three files that contain synchronized segments from an Illumina MiSeq instrument. Pheniqs will assemble an input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input file. Relative input and output file paths are resolved against the working directory which defaults to where you execute pheniqs. You may optionally specify the `base input url` and `base output url` directives.

>```json
{
    "input": [
        "000000000-BDGGG_Lane1_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_R2_001.fastq.gz"
    ]
}
```
>**Example 1.1** Declaring an input read that is [split](glossary.html#split_file_layout) over three gzip compressed FASTQ files.
{: .example}

**Example 1.1** is already a complete and valid Pheniqs configuration! Since no manipulation instructions are specified reads are simply interleaved to the output. Since output is not explicitly declared it defaults to the SAM format and written to standard output.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     FI:i:1  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   13      *       0       0       *       *       0       0       GGACTCCT        B@CCCFC<        FI:i:2  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     FI:i:3  TC:i:3
```
>**Example 1.2** Output header and first 3 records (one complete read) from [interleaving](glossary.html#interleaved_file_layout) the three read segments verbatim into a single SAM formatted stream written to standard output using the configuration file in **Example 1.1**.
{: .example}

# Declaring Output
This simple example is very useful for interleaving raw split read segments into a single CRAM file. CRAM files are the latest indexed and compressed binary encoding of SAM implemented in [HTSlib](glossary.html#htslib) and often provide more efficient compression than the ubiquitous gzip compressed FASTQ format while being much faster to interact with. Packaging your reads in a CRAM container also makes archiving your raw data simple. To write the interleaved output to a compressed CRAM file simply add an `output` directive to **Example 1.1**.

>```json
{
    "input": [
        "000000000-BDGGG_Lane1_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_R2_001.fastq.gz"
    ],
    "output": [ "000000000-BDGGG_raw.cram" ]
}
```
>**Example 1.3** Interleaving three read segments verbatim into a single CRAM file. CRAM files are often much faster to read and write, especially in highly parallelized environments, and also support a rich metadata vocabulary.
{: .example}

# Output Manipulation
The `transform` directive can be used to manipulate the structure of the output read. If omitted all segments of the input are written verbatim to the output, as seen in **Example 1.1** and **Example 1.3**. Since the second segment contains only a technical sequence, and we do not want to write it to the output, we add a `transform` directive to construct an output read from only the first and third segments of the input.

>```json
{
    "input": [
        "000000000-BDGGG_Lane1_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_R2_001.fastq.gz"
    ],
    "transform": { "token": [ "0::", "2::" ] }
}
```
>**Example 1.4** Adding a transform directive composing the output read from only the untouched first and third input segments. Input segments in pheniqs are indexed and referenced using a [zero based coordinate system](glossary.html#zero_based_coordinate) so the first segment is 0.
{: .example}

The [token patterns](manual.html#tokenization) declared in the `token` array of the `transform` directive are made of 3 colon separated integers. The first is the [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment). The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. The two colons are always mandatory. Pheniqs token pattern can address segments from either the 5' (left) or 3' (right) end. Since they mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax they are fairly easy to test.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
```
>**Example 1.5** Output header and first 2 records (one complete read) from interleaving the three read segments to a SAM formatted stream using the configuration file in **Example 1.4**. Notice that only the first and third segments were written to the output.
{: .example}

# Demultiplexing
The reads in our files were sequenced from DNA from 5 individually prepared libraries that were tagged with an 8bp technical sequence before they were pooled together for sequencing. To classify them into read groups we need to examine the first 8 nucleotides of the second segment. Decoding the barcode can be as trivial as comparing two strings if we were absolutely confident no errors occurred during sequencing. However in a real world scenario each nucleotide reported by a sequencing instrument is accompanied by an estimate of the probability the base was incorrectly called. We refer to such an uncertain sequence as an **observed sequence**. The `multiplex` directive is used to declare a single decoder used to classify the reads by examining the segments constructed by the embedded `transform` directive and comparing them to the sequence segments declared in the embedded `codec` directive.

>```json
{
    "multiplex": {
        "transform": { "token": [ "1::8" ] },
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
>**Example 1.6** A `multiplex` decoder directive declaration using the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The embedded `transform` directive is used to extract observed segments from the raw read while the `codec` directive names the possible barcode sequences we expect to find.
{: .example}

In this example we declare a `multiplex` directive that uses the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) algorithm. This algorithm will choose a barcode using a maximum likelihood estimate and reject any classification with a decoding confidence lower than the `confidence threshold` parameter. The `noise` parameter is the prior probability that an observed sequence has not originated from any of the provided barcodes and is often set to the amount of [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution for reads sequenced on the Illumina platform. In the embedded `codec` directive we provide a discrete set of possible decoding results. All `barcode` segment arrays must match the layout declared in the embedded `transform` directive, in this example one 8bp segment. The keys of the `codec` directive can be any unique string. In this example we used the unique barcode nucleotide sequence prefixed with an @ character to remind us this is simply a unique identifier.

>```json
{
    "input": [
        "000000000-BDGGG_Lane1_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_Lane1_S1_L001_R2_001.fastq.gz"
    ],
    "transform": { "token": [ "0::", "2::" ] },
    "multiplex": {
        "transform": { "token": [ "1::8" ] },
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
    },
    "CN": "CGSB",
    "DT": "2018-02-25T07:00:00+00:00",
    "PI": "300",
    "PL": "ILLUMINA",
    "PM": "miseq",
    "SM": "trinidad",
    "flowcell id": "000000000-BDGGG",
    "flowcell lane number": 1
}
```
>**Example 1.7** A complete instruction for demultiplexing with one 8bp barcode segment present on the second input segment to an [interleaved](glossary.html#interleaved_file_layout) SAM stream. The [CN](glossary.html#cn_auxiliary_tag), [DT](glossary.html#dt_auxiliary_tag), [PI](glossary.html#pi_auxiliary_tag), [PL](glossary.html#pl_auxiliary_tag), [PM](glossary.html#pm_auxiliary_tag) and [SM](glossary.html#sm_auxiliary_tag) tags are declared globally and will be added to all read groups, while the [LB](glossary.html#lb_auxiliary_tag) is declared individually for each read group.
{: .example}

The SAM [RG](glossary.html#rg_auxiliary_tag) header tag can contain additional metadata. Capitalized 2 letter directives in a pheniqs configuration file often refer to their corresponding SAM tag and can be either specified for an individual read group or globally for inclusion in all read groups.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:000000000-BDGGG:1:AGGCAGAA     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 5   PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:AGGCAGAA     SM:trinidad
@RG     ID:000000000-BDGGG:1:CGTACTAG     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 4   PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:CGTACTAG     SM:trinidad
@RG     ID:000000000-BDGGG:1:GGACTCCT     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 9   PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:GGACTCCT     SM:trinidad
@RG     ID:000000000-BDGGG:1:TAAGGCGA     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 1   PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:TAAGGCGA     SM:trinidad
@RG     ID:000000000-BDGGG:1:TCCTGAGC     CN:CGSB DT:2018-02-25T07:00:00+00:00    LB:trinidad 8   PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:TCCTGAGC     SM:trinidad
@RG     ID:000000000-BDGGG:1:undetermined CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:000000000-BDGGG:1:undetermined SM:trinidad
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     RG:Z:000000000-BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     RG:Z:000000000-BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
M02455:162:000000000-BDGGG:1:1101:10000:12232   77      *       0       0       *       *       0       0       GTATAGGGGTCACATATAGTTGGTGTGCTTTGTGAACTGCGATCTTGACGG     CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     RG:Z:000000000-BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:CCCCCGGG   XB:f:1.56086e-06
M02455:162:000000000-BDGGG:1:1101:10000:12232   141     *       0       0       *       *       0       0       GTCCTATCCTACTCGGCTTCTCCCCATTTTTCAGACATTTTCCTATCAGTC     CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     RG:Z:000000000-BDGGG:1:GGACTCCT   BC:Z:GGACTCCT   QT:Z:CCCCCGGG   XB:f:1.56086e-06
```
>**Example 1.8** Output header and first 4 records (two complete reads) from demultiplexing using the configuration file in **Example 1.7**. Notice how tags declared globally were added to all read groups. The [RG](glossary.html#rg_auxiliary_tag) and [PU](glossary.html#pu_auxiliary_tag) read group identifiers default to the convention set by [GATK](https://software.broadinstitute.org/gatk/guide/article?id=6472) if you provide the flowcell related directives. The [XB](glossary.html#xb_auxiliary_tag) tag reports the probability the read was incorrectly classified.
{: .example}

# Providing a Prior
If the 5 libraries were pooled in non uniform concentrations we will expect the portion of reads classified to each read group to match those proportions. A prior on the barcode prevalence distribution can be provided for each possible code declared in the `codec` directive using the `concentration` parameter. For convenience the priors do not have to be specified as normalized probabilities. Pheniqs will normalize them when compiling the instructions to sum up to 1.0 minus the value of the `noise` parameter.

>```json
{
    "multiplex": {
        "transform": { "token": [ "1::8" ] },
        "codec": {
            "@AGGCAGAA": { "barcode": [ "AGGCAGAA" ], "concentration": 2 },
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
>**Example 1.6** Adding a prior to one of the code words when declaring a `multiplex` decoder directive. Since the priors are automatically normalized and default to 1, this declaration effectively states that we expect twice as many reads to be classified to @AGGCAGAA than the other 4 read groups.
{: .example}

# Minimum Distance Decoding
Pheniqs can also be instructed to decode the barcodes using the traditional [minimum distance decoder](glossary.html#minimum_distance_decoding), that only consults the edit distance between the expected and observed sequence, by setting the multiplex decoder `algorithm` directive to `mdd`. The MDD decoder however ignores the error probabilities provided by the sequencing instrument and does not compute or report the classification error probability. It is provided for legacy purposes but PAMLD will yield superior results in almost every real world scenario.

# More Efficient Demultiplexing
As we mentioned before reading input from CRAM input can be vastly superior to reading from gzip compressed FASTQ files. If your first step was to package the 3 split FASTQ files into a CRAM file, as shown in **Example 1.3**, you can use that file as input.

>```json
{
    "input": [
        "000000000-BDGGG_raw.cram",
        "000000000-BDGGG_raw.cram",
        "000000000-BDGGG_raw.cram"
    ],
    "transform": { "token": [ "0::", "2::" ] },
    "multiplex": {
        "transform": { "token": [ "1::8" ] },
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
    },
    "CN": "CGSB",
    "DT": "2018-02-25T07:00:00+00:00",
    "PI": "300",
    "PL": "ILLUMINA",
    "PM": "miseq",
    "SM": "trinidad",
    "flowcell id": "000000000-BDGGG",
    "flowcell lane number": 1
}
```
>**Example 1.9** modifying **Example 1.7** to take the CRAM file created in **Example 1.3** as input. When declaring interleaved input you specify the file path as many times as the interleaving resolution. The interleaving resolution is the number of consecutive fragments of the same read that have been interleaved into the file. In this example we expect every read to have 3 consecutive segments.
{: .example}

The output from **Example 1.9** will be identical to the output from **Example 1.7** shown in **Example 1.8** but unlike **Example 1.7** the decoding speed will scale linearly with the number of computational cores available until the system's I/O throughput becomes saturated.
