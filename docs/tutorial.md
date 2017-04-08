<!-- 
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
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
        <li><a                  href="/pheniqs/">Home</a></li>
        <li><a class="active"   href="/pheniqs/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Quick Tutorial
{:.page-title}

* placeholder
{:toc }

The Pheniqs command line interface accepts a [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration file containing a number of separate sections specifying directives for input and output layout, parsing read segments and run parameters. In the [workflow page](workflow.html) you will find some annotated examples of complete configuration files. If you are new to JSON, a [validator](appendix.html#json-validation) can be instrumental for troubleshooting syntax errors. Some parameters are also exposed as command line arguments that override their corresponding configuration file values. A brief description of the command line parameters Pheniqs accepts is always available with the `-h/--help` flags. If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to [install the bundled command line completion](cli.html#zsh-completion) script for a more interactive command line experience.

Pheniqs achieves arbitrary read manipulation in two steps: [tokenization](#tokenization) and [construction](#construction). In the tokenization step you define token patterns that extract a token from an [input segment](glossary.html#input_segment). In the construction step you reference the tokens in [transform patterns](manual.html#transform-pattern) to construct either an [output segment](glossary.html#output_segment), a [multiplex barcode](glossary.html#multiplex_barcode) segment or a [molecular barcode](glossary.html#molecular_barcode).

# File Format

Both input and output [FASTQ](glossary.html#fastq), [SAM, BAM and CRAM](glossary.html#htslib) encoded files are supported. FASTQ can be either uncompressed or gzip compressed. The configuration syntax is sufficiently flexible to easily manipulate [split](glossary.html#split_file_layout), [interleaved](glossary.html#interleaved_file_layout) and [combined](glossary.html#combined_file_layout) file layouts.

# Layout

To declare a layout you define the top level configuration `input`, `token`, `template` and `channel` arrays. If you are demultiplexing you should also define the `multiplex barcode` array. To extract molecular barcodes you can use the `molecular barcode` property.

## Input

The top level configuration `input` array is an ordered list of input file paths that Pheniqs can use as input feeds. Pheniqs assembles an input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input feed.

>```json
{
    "input": [
        "HK5NHBGXX_l01n01.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n03.fastq.gz", 
        "HK5NHBGXX_l01n04.fastq.gz"
    ]
}
```
>Assembling an input read from four segments [split](glossary.html#split_file_layout) over four gzip compressed FASTQ files.
{: .example}


[Interleaved](glossary.html#interleaved_file_layout) files contain multiple consecutive segments of the same read. To assemble an input read from an interleaved feed simply repeat the path to reference the same feed multiple times, once for each [segment](glossary.html#segment).

>```json
{
    "input": [
        "HK5NHBGXX_l01.bam", 
        "HK5NHBGXX_l01.bam", 
        "HK5NHBGXX_l01.bam", 
        "HK5NHBGXX_l01.bam"
    ]
}
```
>Constructing an input read from four segments [interleaved](glossary.html#interleaved_file_layout) into a single BAM file.
{: .example}

## Tokenization

A [token pattern](manual.html#tokenization) is made of 3 colon separated integers. The first is the [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment). The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment.

>```json
{
    "token": [
        "0::", 
        "1::8", 
        "2::8", 
        "3::"
    ]
}
```
>**Tokenizing a standard dual-indexed paired-end Illumina run.** Reads have two [template](glossary.html#template_sequence) sequences in the first and fourth segments and two 8bp [multiplex barcodes](glossary.html#multiplex_barcode) in the second and third segments. Segments are often sequenced one cycle more than required to ensure good quality on the last relevant cycle and so although the second and third sequences contain 9 cycles, only the first 8 are the desired [technical sequence](glossary.html#template_sequence). In this example we preserve the complete template sequences, which can be of a variable length.
{: .example}

## Construction

A [transform pattern](manual.html#transform-pattern) is made of one or more token references separated by the **:** concatenation operator. To construct an [output segment](glossary.html#output_segment) declare a transform pattern in the top level configuration `template` array. To construct a [multiplex](glossary.html#multiplex_barcode) or [molecular](glossary.html#molecular_barcode) barcode segment declare a a transform pattern in the top level configuration `multiplex barcode` or `molecular barcode` array, respectively.

>```json
{
    "token": [
        "0:6:", 
        "1::8", 
        "2::8", 
        "3::",
        "0::6"
    ],
    "template": [
        "0", 
        "3"
    ],
    "multiplex barcode": [
        "1", 
        "2"
    ],
    "molecular barcode": [
        "4"
    ]
}
```
Extending the previous example, we further assume the first 6 cycles of the first input segment contain a single molecular barcode. Output template reads have two segments constructed from tokens 0 and 3. Tokens 1 and 2 are used to construct the two [technical sequences](glossary.html#template_sequence) of the multiplex barcode. A molecular barcode will be constructed from token 4. Notice that token 0 only begins with the 7th cycle so in this case the molecular barcode technical sequence will not be present on the first segment of the output read.
{: .example}

## Output

To keep things simple we assume for now that every [library](glossary.html#library) in your experiment maps to a different [read group](glossary.html#read_group), but you will later see that Pheniqs can support even more flexible layouts. For every library you declare a channel specification element in the top level configuration `channel` array. A channel specification element defines read group classification, output layout and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). Multiple channels can write to the same concrete output file to create a [combined](glossary.html#combined_file_layout) output layout. Every read is guaranteed to have all its segments written contiguously and in-order.

The nested channel specification `barcode` array is an ordered list of multiplex barcode [technical sequences](glossary.html#technical_sequence). Degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) bases are currently not allowed in a multiplex barcode sequence and will be replaced with the **N** (any nucleotide) symbol before decoding. The size of the `barcode` array must be the same as the size to the top level `multiplex barcode` array and corresponds to the number of multiplex barcode sets.

The nested channel specification `output` array is an ordered list of output file paths that Pheniqs can use to as output feeds. Pheniqs writes one output read [segment](glossary.html#segment) to each output feed. When writing an [interleaved](glossary.html#interleaved_file_layout) output simply repeat the path to reference the same feed multiple times, once for each [segment](glossary.html#segment). The size of the `output` array must be the same as the size of the top level `template` array and corresponds to the number of segments in an output read. One syntactical shortcut to this rule is when all segments are interleaved into the same feed, in which case you may specify the path only once.

The nested channel specification `concentration` property can be used to provide a prior probability of observing a read from this channel, most commonly the relative pooling concentration of the multiplexed library.

Each channel specification defines a [read group](glossary.html#read_group). A corresponding [@RG header tag](glossary.html#rg_header_tag) will be added to the header of every output feed. Reads in the read group will be tagged with the corresponding [RG auxiliary tag](glossary.html#rg_auxiliary_tag).

If you wish to retain undetermined reads you must explicitly define an undetermined channel and set the nested channel specification `undetermined` property to `true`. You may also set some header properties and auxiliary tags [globally](manual.html#global-auxiliary-tags) and they will be assigned to all channels.

>```json
{
    "CN": "NYU CGSB", 
    "DT": "2016-07-15T07:00:00+00:00",
    "PI": "1500",
    "PL": "ILLUMINA",
    "PM": "HiSeq 2500",
    "channel": [
        {
            "RG": "HCJFNBCXX:1:TAAGGCGATCTTACGC", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "DS": "TAAGGCGATCTTACGC read group description", 
            "LB": "TAAGGCGATCTTACGC read group library name", 
            "SM": "TAAGGCGATCTTACGC read group sample name", 
            "barcode": [
                "TAAGGCGA", 
                "CTCTCTAT"
            ], 
            "output": [
                "HG7CVAFXX_l01s01_TAAGGCGATCTTACGC.cram"
            ],
            "concentration": 1
        },
        {
            "RG": "HCJFNBCXX:1:undetermined", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "output": [
                "HG7CVAFXX_l01s01_undetermined.fastq.gz", 
                "HG7CVAFXX_l01s02_undetermined.fastq.gz"
            ], 
            "undetermined": true
        }
    ]
}
```
> Defining one standard output channel and one for undetermined reads. Notice that although all channels expect a read with the same number of segments, each can write it to a different layout. In this example the first channel writes to an [interleaved](glossary.html#interleaved_file_layout) CRAM file while the second [splits](glossary.html#split_file_layout) the read to two FASTQ files. The global [CN](glossary.html#cn_auxiliary_tag), [DT](glossary.html#dt_auxiliary_tag), [PI](glossary.html#pi_auxiliary_tag), [PL](glossary.html#pl_auxiliary_tag) and [PM](glossary.html#pm_auxiliary_tag) properties will be added to all read groups.
{: .example}

# Decoding

Pheniqs implements a novel [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that directly estimates the decoding likelihood from the base calling error probabilities provided by the sequencing platform. When output is SAM encoded Pheniqs will write the multiplex barcode decoding error probability to the [DQ](glossary.html#dq_auxiliary_tag) auxiliary tag. Pheniqs also implements a traditional [minimum distance decoder](glossary.html#minimum_distance_decoding) that consults the edit distance between the expected and observed sequence. To select a decoder set the top level configuration `decoder` property to either `pamld` or `mdd`. The default decoder is `pamld`.

>```json
{ "decoder": "mdd" }
```
> The default decoder is the phred-adjusted maximum likelihood decoder. Set the `decoder` property to `mdd` to use minimum distance decoding.
{: .example}

When using phred-adjusted maximum likelihood decoding you may want to set the `confidence` and `noise` properties. The `confidence` property is the threshold on the decoding probability for the decoder to declare a successful classification. If the decoder fails to classify the read it is considered undetermined. The value of `confidence` defaults to **0.99** but in practice depends on the application. Since the multiplex barcode decoding error probability is written to the [DQ](glossary.html#dq_auxiliary_tag) auxiliary tag you can set `confidence` to a relatively lax value and leave the decision to exclude the read for downstream analysis.

>```json
{ "confidence": 0.99 }
```
> Setting the confidence to 0.99 is equivalent to a 1% error probability.
{: .example}

The `noise` property can be used to provide a prior probability of observing an undetermined read. For the Illumina platform the most common value will be the concentration of the [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution, which is usually 1% and corresponds to a value of **0.01**. The `noise` property defaults to **0**.

>```json
{ "noise": 0.01 }
```
> Setting the noise prior to 0.01 when 1% of PhiX control was spiked in to the solution.
{: .example}

# Performance Tuning

Several parameters can help [improve performance](manual.html#performance-tuning) on multi core machines. Setting the `threads` property to the number of available cores will use a benchmark based heuristic algorithm to set the more low level properties and should be a good starting point.