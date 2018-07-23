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
>Output from interleaving the three read segments into a single unmanipulated SAM formatted stream written to standard output.
{: .example}

This simple configuration can be very useful for repacking split read segments into a single CRAM file, which are much faster to read from. Packing your reads into a single file also makes to archiving easier and should usually be your first step. To write the interleaved output to a compressed CRAM file simply add an output directive.

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
>Interleaving three read segments into a single efficient CRAM file without any manipulation. CRAM files are often much faster to read and write, especially in highly parallelized environments, and also offer a rich metadata vocabulary.
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
>Adding a template directive composing an output read from the untouched first and third input segments. Input segments in pheniqs are indexed and referenced using a [zero based coordinate system](glossary.html#zero_based_coordinate) so the first segment is 0.

The [token patterns](manual.html#tokenization) declared in the token array of the template directive are made of 3 colon separated integers. The first is the [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment). The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. The two colons are always mandatory. Pheniqs token pattern can address segments from either the 5' (left) or 3' (right) end and mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
```
>>Output from interleaving the three read segments into a SAM formatted stream. Notice that only the first and third segments were written.
{: .example}

## Read Group classification
The reads in our files are DNA sequenced from 6 different libraries that were tagged with an 8bp technical sequence and pooled together for sequencing. We want to classify the reads into read groups by examining the first 8 nucleotides of the second segment. Decoding the barcode would have been a simple matter of comparing two words if we were absolutely confident no errors were made during sequencing, which we know is not true. The top level `multiplex` element declares the multiplex decoder.

>```json
{
    "multiplex": {
        "algorithm": "pamld",
        "noise": 0.02,
        "confidence threshold": 0.95,
        "template": { "token": [ "1::8" ] },
        "codec": {
            "@AGGCAGAA": { "barcode": [ "AGGCAGAA" ] },
            "@CGTACTAG": { "barcode": [ "CGTACTAG" ] },
            "@GGACTCCT": { "barcode": [ "GGACTCCT" ] },
            "@TAAGGCGA": { "barcode": [ "TAAGGCGA" ] },
            "@TCCTGAGC": { "barcode": [ "TCCTGAGC" ] }
        }
    }
}
```

In this example we declare the decoding algorithm and the parameters it needs: an expected noise level (portion of reads that should not be classified to any of the barcodes) and confidence threshold. We also provide the already familiar template element, in this context to instruct Pheniqs what portion of the input read we want to decode a barcode from. In the `codec` element we name the possible values we expect to find. The `barcode` element in each of the options must match in layout to what you declared in the decoder's `template`. Degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) bases are currently not allowed in a multiplex barcode sequence.

>```json
{
    "CN": "CGSB",
    "DT": "2018-02-25T07:00:00+00:00",
    "PI": "300",
    "PL": "ILLUMINA",
    "PM": "miseq",
    "flowcell id": "BDGGG",
    "input": [
        "BDGGG_s01.fastq.gz",
        "BDGGG_s02.fastq.gz",
        "BDGGG_s03.fastq.gz"
    ],
    "output": [ "BDGGG.cram" ],
    "template": { "token": [ "0::", "2::" ] },
    "multiplex": {
        "algorithm": "pamld",
        "noise": 0.02,
        "confidence threshold": 0.95,
        "template": { "token": [ "1::8" ] },
        "codec": {
            "@AGGCAGAA": { "barcode": [ "AGGCAGAA" ] },
            "@CGTACTAG": { "barcode": [ "CGTACTAG" ] },
            "@GGACTCCT": { "barcode": [ "GGACTCCT" ] },
            "@TAAGGCGA": { "barcode": [ "TAAGGCGA" ] },
            "@TCCTGAGC": { "barcode": [ "TCCTGAGC" ] }
        }
    }
}
```
> Demultiplexing with one 8bp barcode on the second segment to an [interleaved](glossary.html#interleaved_file_layout) CRAM file. The [CN](glossary.html#cn_auxiliary_tag), [DT](glossary.html#dt_auxiliary_tag), [PI](glossary.html#pi_auxiliary_tag), [PL](glossary.html#pl_auxiliary_tag) and [PM](glossary.html#pm_auxiliary_tag) tags are declared globally and will be added to all read groups.
{: .example}

Additional read group metadata can be declared in the individual codec elements or be placed in the top level configuration element to be included in all read groups. Output from executing the above configuration will look something like this

```
@HD     VN:1.0  SO:unknown      GO:query
@PG     ID:pheniqs      PN:pheniqs      CL:pheniqs demux --config test/BDGGG/BDGGG_interleave.json      VN:2.0.3-beta-7-g3bb189989b86be01c266719e9b449367cf4157e5
@RG     ID:BDGGG:AGGCAGAA       CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:AGGCAGAA
@RG     ID:BDGGG:CGTACTAG       CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:CGTACTAG
@RG     ID:BDGGG:GGACTCCT       CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:GGACTCCT
@RG     ID:BDGGG:TAAGGCGA       CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:TAAGGCGA
@RG     ID:BDGGG:TCCTGAGC       CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:TCCTGAGC
@RG     ID:BDGGG:undetermined   CN:CGSB DT:2018-02-25T07:00:00+00:00    PI:300  PL:ILLUMINA     PM:miseq        PU:BDGGG:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     RG:Z:BDGGG:GGACTCCT     BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     RG:Z:BDGGG:GGACTCCT     BC:Z:GGACTCCT   QT:Z:B@CCCFC<   XB:f:1.56496e-06
```

# Decoding

Pheniqs implements a novel [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that directly estimates the decoding likelihood from the base calling error probabilities provided by the sequencing platform. When output is SAM encoded Pheniqs will write the multiplex barcode decoding error probability to the [XB](glossary.html#dq_auxiliary_tag) auxiliary tag. Pheniqs also implements a traditional [minimum distance decoder](glossary.html#minimum_distance_decoding) that consults the edit distance between the expected and observed sequence. To select a decoding algorithm set `algorithm` property to either `pamld` or `mdd`.

When using phred-adjusted maximum likelihood decoding you may want to set the `confidence threshold` and `noise` properties. The `confidence threshold` property is the threshold on the decoding probability for the decoder to declare a successful classification. If the decoder fails to classify the read it is considered undetermined. The value of `confidence threshold` defaults to **0.99** but in practice depends on the application. Since the multiplex barcode decoding error probability is written to the [XB](glossary.html#dq_auxiliary_tag) auxiliary tag you can set `confidence threshold` to a relatively lax value and leave the decision to exclude the read for downstream analysis.

The `noise` property can be used to provide a prior probability of observing an undetermined read. For the Illumina platform the most common value will be the concentration of the [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution, which is usually 1% and corresponds to a value of **0.01**. The `noise` property defaults to **0.01**.
