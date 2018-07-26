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
        <li><a                  href="/pheniqs/1.0/">Home</a></li>
        <li><a                  href="/pheniqs/1.0/tutorial.html">Tutorial</a></li>
        <li><a class="active"   href="/pheniqs/1.0/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/1.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/1.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/1.0/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/1.0/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Pheniqs Documentation
{:.page-title}

* placeholder
{:toc}

The Pheniqs command line interface accepts a [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration file. Some parameters are also exposed as command line arguments that override their corresponding configuration file values. The configuration file contains a number of separate sections specifying directives for input and output layout, parsing read segments and run parameters. The different sections of the configuration file are described bellow. The [workflow page](workflow.html) contains some annotated examples of complete configuration files. 

Pheniqs achieves arbitrary read manipulation in two steps: [tokenization](#tokenization) and [construction](#construction). In the tokenization step you define token patterns that extract a token from an [input segment](glossary.html#input_segment). In the construction step you reference the tokens in [transform patterns](#transform-pattern) to construct either an [output](glossary.html#output_segment), a [multiplex barcode](glossary.html#multiplex_barcode) or a [molecular barcode](glossary.html#molecular_barcode) segment.

# File Format
Both input and output [FASTQ](glossary.html#fastq), [SAM, BAM and CRAM](glossary.html#htslib) encoded files are supported. FASTQ can be either uncompressed or gzip compressed. The configuration syntax is sufficiently flexible to manipulate [split](glossary.html#split_file_layout), [interleaved](glossary.html#interleaved_file_layout) and [combined](glossary.html#combined_file_layout) file layouts.

# Layout

To declare a layout you define the top level configuration `input`, `token`, `template` and `channel` arrays. If you are demultiplexing you should also define the `multiplex barcode` array. To extract molecular barcodes you can use the `molecular barcode` array.


## Input

The top level configuration `input` array is an ordered list of input file paths. Pheniqs assembles an input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input feed.

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
>
>```json
{
    "input": [
        "HK5NHBGXX_l01n01.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n01.fastq.gz" 
    ]
}
```
> Constructing an input read from four segments. Segments 0 and 3 are [interleaved](glossary.html#interleaved_file_layout) in `HK5NHBGXX_l01n01.fastq.gz` and segments 1 and 2 are interleaved into `HK5NHBGXX_l01n02.fastq.gz`.
{: .example}

## Tokenization

A token pattern is made of 3 colon separated integers. The first is the mandatory [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment) enumerated by `input`. The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. **start** coordinate and **end** coordinate can take positive or negative values to access the segment from either the 5' (left) or 3' (right) end and mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax. The two colons are always mandatory.

>```json
{
    "token": [
        "0:0:", 
        "1:0:8"
    ]
}
```
>
>The first token spans the entire first input segment. The second spans the first 8 cycles of the second segment.
{: .example}

>| **Cycle**      | `012345678` |
>| **Nucleotide** | `GGACTCCTA` |
>
>| Pattern   | Segment   | Start   | End   | Cycles      | Bitmap      | Description                           |
>| :-------- | :-------- | :------ | :---- | :---------- | :---------- | :------------------------------------ |
>| `0:0:3`   | `0`       | `0`     | `3`   | `GGA`       | `+++------` | First 3 cycles of segment 0           |
>| `0:0:`    | `0`       | `0`     | `9`   | `GGACTCCTA` | `+++++++++` | All of segment 0                      |
>| `0::`     | `0`       | `0`     | `9`   | `GGACTCCTA` | `+++++++++` | All of segment 0                      |
>| `0:-4:`   | `0`       | `5`     | `9`   | `CCTA`      | `-----++++` | Last 4 cycles of segment 0            |
>| `0:-5:-1` | `0`       | `4`     | `8`   | `TCCT`      | `----++++-` | Cycles 1 to 5 from the 3' end         |
>| `0:3:-2`  | `0`       | `3`     | `7`   | `CTCC`      | `---++++--` | All but the first 3 cycles and last 2 |
>
> Some examples of tokenization of a hypothetical segment with index 0 and 9 cycles.
{: .example}

## Transform pattern

A transform pattern is made of one or more token references separated by the **:** concatenation operator. A token reference is the [zero based](glossary.html#zero_based_coordinate) index of the token pattern in the `token` array. Appending the left hand side reverse complementarity **~** operator will concatenate the reverse complemented sequence of the token. Each token reference is evaluated before concatenation so **~** evaluation precedes **:** evaluation.

>| Pattern | Description                                                                                        |
>| ------- | :------------------------------------------------------------------------------------------------- |
>| `0`     | The simplest possible transform will assemble a segment from just token 0.                         |
>| `~0`    | Assemble a segment from the reverse complemented token 0.                                          |
>| `0:1`   | Assemble a segment by concatenating token 0 and token 1                                            |
>| `~0:1`  | Assemble a segment by concatenating the reverse complement of token 0 and token 1.                 |
>
> Several examples of transform patterns for constructing new segments. 
{: .example}

Notice that a negative token **start** coordinate is equivalent to the corresponding positive token **end** coordinate on the reverse complemented strand and vice verse. More formally the token `0:-x:-y` is equivalent to `0:y:x` if applied to the reverse complement of the segment. For instance to concatenate to the first output segment the first 6 bases of the reverse complemented strand of the first input segment you would define token `0:-6:` and then reference it in the transform pattern for output segment 0 as `~0`.

## Construction

To construct an [output segment](glossary.html#output_segment) declare a transform pattern in the top level configuration `template` array. To construct a [multiplex barcode](glossary.html#multiplex_barcode) segment declare a a transform pattern in the top level configuration `multiplex barcode` array. To construct a [molecular barcode](glossary.html#molecular_barcode) segment declare a a transform pattern in the top level configuration `molecular barcode` array.

The size of the `template` and `multiplex barcode` arrays implicitly defines the number of segments in an output read and the number of barcode sets, respectively.

>```json
{
    "token": [
        "0:6:", 
        "3::-6",
        "1::8", 
        "2::8",
        "0::6",
        "3:-6:"
    ],
    "template": [
        "0", 
        "1"
    ],
    "multiplex barcode": [
        "2", 
        "3"
    ],
    "molecular barcode": [
        "4",
        "~5"
    ]
}
```
>Consider a dual indexed paired end read with two 8 cycle multiplex barcode segments on the second and third input segments and two molecular barcodes, with one on the first 6 cycles of the first segment and another, reverse complemented, on the last 6 cycles of the last segment. We want to extract both molecular barcodes and remove the them from the segments.
{: .example}

## Output

Output layout is defined by the top level `channel` configuration array. Each channel element defines [read group](glossary.html#read_group) classification, output layout and [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf). If you wish to retain undetermined reads you must explicitly define an undetermined channel and set the nested channel specification `undetermined` property to `true`. You may also set some of the read group properties [globally](#global-auxiliary-tags) and they will be assigned to all read group decelerations that do not already provide a value for those properties.

### Read Group

Read groups can be either explicitly declared in the top level `read group` array or implicitly in the channel element. Separating the read group declaration from the channel is beneficial when multiple channels are assigned the same read group. The `read group` array is interpreted before the `channel` array and both are interpreted in the order the elements are present in the configuration file. If a channel `RG` property does not reference a previously declared read group, either implicitly or explicitly, it is assumed to be an implicit read group declaration and read group properties are decoded from the channel element, so both read group declaration styles can coexist. 

Explicitly declared [read groups](glossary.html#read_group) are defined in the the top level `read group` configuration array. A channel can reference the read group by setting its nested `RG` property to the value of the `ID` property. 
 

| Name                                      | Description                                      | Type   |
| :---------------------------------------- | :----------------------------------------------- | :----- |
| **ID**                                    | Read group identifier                            | string |
| **[LB](glossary.html#lb_auxiliary_tag)**  | Library name                                     | string |
| **[SM](glossary.html#sm_auxiliary_tag)**  | Sample name                                      | string |
| **[PU](glossary.html#pu_auxiliary_tag)**  | Platform unit unique identifier                  | string |
| **CN**                                    | Name of sequencing center producing the read     | string |
| **DS**                                    | Description                                      | string |
| **DT**                                    | ISO8601 date the run was produced                | string |
| **PI**                                    | Predicted median insert size                     | string |
| **[PL](glossary.html#pl_auxiliary_tag)**  | Platform or technology used to produce the reads | string |
| **PM**                                    | Platform model                                   | string |
| **PG**                                    | Programs used for processing the read group      | string |


The **PL** property is one of *CAPILLARY*, *LS454*, *ILLUMINA*, *SOLID*, *HELICOS*, *IONTORRENT*, *ONT*, *PACBIO*

The **ID** is mandatory and must be unique among all read groups.

>```json
{
    "read group": [
        {
            "ID": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "LB": "TAAGGCGATCTTACGC read group library name", 
            "SM": "TAAGGCGATCTTACGC read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "TAAGGCGATCTTACGC read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500"
        },
        {
            "ID": "HK5NHBGXX:1:undetermined", 
            "LB": "undetermined read group library name", 
            "SM": "undetermined read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "undetermined read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500"
        }
    ]
}
```
> Declaring one read group for a multiplexed channel and one for undetermined reads.
{: .example}

### Channel
Multiple channels can write to the same concrete output file to create a [combined](glossary.html#combined_file_layout) output layout. Channels are guaranteed to write all segments of a read contiguously and in-order. The order of the reads in the output is **not** guaranteed to be consistent with their order in the input or to be stable between successive executions.

The nested channel `barcode` array is an ordered list of multiplex barcode [technical sequences](glossary.html#technical_sequence). Degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) bases are currently not allowed in the [multiplex barcode](glossary.html#multiplex_barcode) sequence and will be replaced with the **N** (any nucleotide) symbol before decoding. The size of the `barcode` array must be the same as the size to the top level `multiplex barcode` array and corresponds to the number of multiplex barcode sets.

The nested channel `output` array is an ordered list of output file paths that Pheniqs can use as output feeds. Pheniqs writes one [segment](glossary.html#segment) from the output read to each output feed. When writing an [interleaved](glossary.html#interleaved_file_layout) output simply repeat the path to reference the same feed multiple times, once for each [segment](glossary.html#segment). The size of the `output` array must be the same as the size of the top level `template` array and corresponds to the number of segments in an output read. One syntactical shortcut to this rule is when all segments are interleaved into the same feed, in which case you may specify the path only once.

The nested channel `concentration` property is only used by the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) and can be used to provide a prior probability of observing a read from this channel. `concentration` values for all channels are normalized so that their normalized sum equals `1.0 - noise`. This allows the user to directly specify the pooled concentration and it will be converted to a corresponding prior probability. The `concentration` element must be either specified on all non-undetermined channels or on none of them. If `concentration` values are not specified they default to **1** on all non undetermined channels and so result in a uniform prior distribution. Notice that unlike `concentration` the `noise` value is specified as a probability value between 0 and 1.

The nested channel `RG` property is a reference to a read group `ID`. A corresponding [@RG header tag](glossary.html#rg_header_tag) will be added to the header of every output feed. Reads in the read group will be tagged with the corresponding [RG auxiliary tag](glossary.html#rg_auxiliary_tag).

Setting the nested channel `undetermined` property to `true` identifies that channel as the destination for reads that fail to be classified. Setting the `undetermined` property to `true` in more than one channel will result in a validation error. The undetermined channel may optionally omit the `barcode` element. If omitted the `undetermined` property defaults to `false`.

>```json
{
    "read group": [
        {
            "ID": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "LB": "TAAGGCGATCTTACGC read group library name", 
            "SM": "TAAGGCGATCTTACGC read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "TAAGGCGATCTTACGC read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500"
        },
        {
            "ID": "HK5NHBGXX:1:undetermined", 
            "LB": "undetermined read group library name", 
            "SM": "undetermined read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "undetermined read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500"
        }
    ],
    "channel": [
        {
            "RG": "HCJFNBCXX:1:TAAGGCGATCTTACGC",
            "barcode": [
                "TAAGGCGA",
                "CTCTCTAT"
            ],
            "concentration": 1,
            "output": [
                "HG7CVAFXX_l01s01_TAAGGCGATCTTACGC.cram"
            ]
        },
        {
            "RG": "HCJFNBCXX:1:undetermined",
            "output": [
                "HG7CVAFXX_l01s01_undetermined.fastq.gz",
                "HG7CVAFXX_l01s02_undetermined.fastq.gz",
            ],
            "undetermined": true
        }
    ]
}
```
> Defining one standard output channel and one for undetermined reads with explicitly declared read groups. Notice how the first channel writes to an [interleaved](glossary.html#interleaved_file_layout) CRAM file while the second [splits](glossary.html#split_file_layout) the read to two FASTQ files.
>
>```json
{
    "channel": [
        {
            "RG": "HCJFNBCXX:1:TAAGGCGATCTTACGC",
            "LB": "TAAGGCGATCTTACGC read group library name", 
            "SM": "TAAGGCGATCTTACGC read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "TAAGGCGATCTTACGC read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500",
            "barcode": [
                "TAAGGCGA",
                "CTCTCTAT"
            ],
            "concentration": 1,
            "output": [
                "HG7CVAFXX_l01s01_TAAGGCGATCTTACGC.cram"
            ]
        },
        {
            "RG": "HCJFNBCXX:1:undetermined",
            "LB": "undetermined read group library name", 
            "SM": "undetermined read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "CN": "NYU CGSB", 
            "DS": "undetermined read group description", 
            "DT": "2016-07-15T07:00:00+00:00",
            "PI": "1500",
            "PL": "ILLUMINA",
            "PM": "HiSeq 2500",
            "output": [
                "HG7CVAFXX_l01s01_undetermined.fastq.gz",
                "HG7CVAFXX_l01s02_undetermined.fastq.gz",
            ],
            "undetermined": true
        }
    ]
}
```
> This declaration is equivalent to the previous one but the read groups are implicitly declared in the channels.
>
```json
{
    "CN": "NYU CGSB", 
    "DT": "2016-07-15T07:00:00+00:00",
    "PI": "1500",
    "PL": "ILLUMINA",
    "PM": "HiSeq 2500",
    "channel": [
        {
            "RG": "HCJFNBCXX:1:TAAGGCGATCTTACGC",
            "LB": "TAAGGCGATCTTACGC read group library name", 
            "SM": "TAAGGCGATCTTACGC read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "DS": "TAAGGCGATCTTACGC read group description", 
            "barcode": [
                "TAAGGCGA",
                "CTCTCTAT"
            ],
            "concentration": 1,
            "output": [
                "HG7CVAFXX_l01s01_TAAGGCGATCTTACGC.cram"
            ]
        },
        {
            "RG": "HCJFNBCXX:1:undetermined",
            "LB": "undetermined read group library name", 
            "SM": "undetermined read group sample name", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTTACGC", 
            "DS": "undetermined read group description", 
            "output": [
                "HG7CVAFXX_l01s01_undetermined.fastq.gz",
                "HG7CVAFXX_l01s02_undetermined.fastq.gz",
            ],
            "undetermined": true
        }
    ]
}
```
> In this equivalent declaration the [CN](glossary.html#cn_auxiliary_tag), [DT](glossary.html#dt_auxiliary_tag), [PI](glossary.html#pi_auxiliary_tag), [PL](glossary.html#pl_auxiliary_tag) and [PM](glossary.html#pm_auxiliary_tag) properties are declared globally and will be added to both read groups.
{: .example}

# Decoding

Pheniqs implements a novel [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that directly estimates the decoding likelihood from the base calling error probabilities provided by the sequencing platform. When output is SAM encoded Pheniqs will write the multiplex barcode decoding error probability to the [DQ](glossary.html#dq_auxiliary_tag) auxiliary tag. Pheniqs also implements a traditional [minimum distance decoder](glossary.html#minimum_distance_decoding) that consults the edit distance between the expected and observed sequence. To select a decoder set the top level configuration `decoder` property to either `pamld` or `mdd`. The default decoder is `pamld`.

>```json
{ "decoder": "mdd" }
```
> The default decoder is the phred-adjusted maximum likelihood decoder. Set the `decoder` property to `mdd` to use minimum distance decoding.
{: .example}

A third "meta" decoder called `benchmark` will decode using both `pamld` and `mdd` and report the results of the comparison in some custom user space SAM auxiliary tags; this option is primarily used for evaluation.

When using the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) you may want to set the `confidence` and `noise` properties. The `confidence` property is the threshold on the decoding probability for the decoder to declare a successful classification. If the decoder fails to classify the read it is considered undetermined. The value of `confidence` defaults to **0.99** but in practice depends on the application. Since the multiplex barcode decoding error probability is written to the [DQ](glossary.html#dq_auxiliary_tag) auxiliary tag you can set `confidence` to a relatively lax value and leave the decision to exclude the read for downstream analysis.

>```json
{ "confidence": 0.99 }
```
> Setting the confidence to 0.99 is equivalent to a 1% error probability.
{: .example}

The `noise` property is only used by the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) and takes a floating point value between 0 and 1. It is the probability that a read is background noise and does not belong to any of the multiplexed libraries. For the Illumina platform the most common value will be the amount of [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution, which is usually 1% and corresponds to a value of **0.01**. Assays with low base diversity sometimes use a higher concentration to compensate and some applications can go as high as 50%. If you know a priori that some amount of sequences present in the pool represent contamination or multiplexed libraries you are not declaring in the configuration, their concentration can be factored into this parameter. The `noise` property defaults to **0**.

>```json
{ "noise": 0.01 }
```
> Setting the noise prior to 0.01 when 1% of PhiX control was spiked in to the solution.
{: .example}

The `distance tolerance` array is a list of integers considered only by the [minimum distance decoder](glossary.html#minimum_distance_decoding). The size of the `distance tolerance` array must be the same as the size of the `multiplex barcode` array and is the maximum edit distance allowed for each barcode set to still be considered a match. If omitted, it defaults to the maximum number of correctable errors computed from the pairwise Hamming distance for each barcode set. If set to a value larger than the maximum number of correctable errors it will be ignored.

>```json
{
    "distance tolerance": [
        2,
        1
    ]
}
```
>Setting a maximum edit distance of 2 for the first barcode and 1 for the second.
{: .example}

The `decoder masking threshold` is considered by both the [minimum distance decoder](glossary.html#minimum_distance_decoding) and the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) but is more suitable for using with minimum distance decoding. If set to a value larger than 0, any cycle on a multiplex barcode with a quality lower than `decoder masking threshold` will be set to **N** before decoding.

>```json
{ "decoder masking threshold": 8 }
```
>Any cycle with a quality of less than 8 will be set to **N**, the "any nucleotide" symbol.
{: .example}

# Relative paths

Setting global path prefixes makes the configuration files more portable. If specified, the `base input path` and `base output path` are used as a prefix to **relative** paths defined in the `input` and `output` elements respectively. A path is considered relative if it **does not** begin with a **/** character.

>```json
{
    "base input path": "/volume/alpha/HK5NHBGXX/lane_01", 
    "input": [
        "HK5NHBGXX_l01n01.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n01.fastq.gz" 
    ]
}
```
>
> will resolve to 
>
>```json
{
    "input": [
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_l01n01.fastq.gz", 
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_l01n02.fastq.gz", 
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_l01n02.fastq.gz", 
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_l01n01.fastq.gz"
    ]
}
```
{: .example}

All paths specified in a configuration file can be made relative to your home directory

>```json
{
    "base input path": "~/HK5NHBGXX", 
    "input": [
        "HK5NHBGXX_l01n01.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n02.fastq.gz", 
        "HK5NHBGXX_l01n01.fastq.gz" 
    ]
}
```
>
> will resolve to 
>
>```json
{
    "input": [
        "/home/lg/HK5NHBGXX/HK5NHBGXX_l01n01.fastq.gz", 
        "/home/lg/HK5NHBGXX/HK5NHBGXX_l01n02.fastq.gz", 
        "/home/lg/HK5NHBGXX/HK5NHBGXX_l01n02.fastq.gz", 
        "/home/lg/HK5NHBGXX/HK5NHBGXX_l01n01.fastq.gz"
    ]
}
```
> The `~` character in the beginning of a path is interpreted as the `HOME` environment variable on most POSIX shells and resolves to the currently logged in user. Setting a base path that begins with `~` makes the configuration file portable between systems where you have a different login.
{: .example}

# Phred offset

The `input phred offset` and `output phred offset` are applicable only to [FASTQ](glossary.html#fastq) files and specify the Phred scale decoding and encoding [offset](https://en.wikipedia.org/wiki/FASTQ_format#Encoding), respectively. The default value for both is **33**, knowns as the [Sanger format](glossary.html#sanger_format). The binary [HTSlib](glossary.html#htslib) formats BAM and CRAM encode the quality value numerically and so require no further manipulation. The [sequence alignment map format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) states that text encoded SAM records always use the Sanger format.

>```json
{
    "output phred offset": 33, 
    "input phred offset": 33
}
```
> Setting the `input phred offset` and `output phred offset` to 33 which is the default.
{: .example}

# Leading Segment

Pheniqs constructs new segments for the output read and must generate corresponding identifiers and replicate metadata from the input to the output segments. Since segments can potentially disagree on metadata, one input segment is elected as the leader and used as a template when constructing the output segments. The leading segment property is a reference to an [input segment index](glossary.html#input_segment) and defaults to **0** if omitted.

>```json
{ "leading segment": 0 }
```
>The leading segment defaults to be the first input segment.
{: .example}

# Pass filter reads

Some reads are marked by the sequencing platform as not passing the vendor quality control. For instance Illumina sequencers perform an internal [quality filtering procedure](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-percent-pf-technical-note-770-2014-043.pdf) called chastity filter, and reads that pass this filter are called PF for pass-filter. This can be [signaled](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) on the comment portion of the read identifier in [FASTQ](glossary.html#fastq) files or the **512** flag on an [HTSlib](glossary.html#htslib) flag. Setting this parameter to **true** will allow Pheniqs to decode those unfiltered reads and include them in the output. If the paramter is set to **false** those reads are dropped. If omitted `include filtered` defaults to **false**.

>```json
{ "include filtered": false }
```
{: .example}

# Global auxiliary tags

Some of the [read group](glossary.html#read_group) header properties and SAM auxiliary tags are often identical in most, if not all, of the read groups you define. To streamline your configuration file you may optionally specify those properties globally in the root element and they will be replicated in every read group that does not explicitly define the property. To clarify, a property defined inside a `read group` or `channel` element takes precedence over the globally defined tags. You may think of a global property as the default value.

The following tags are allowed globally in the root element:

| Name     | Description                                      | Type   |
| :------- | :----------------------------------------------- | :----- |
| **CN**   | Name of sequencing center producing the read     | string |
| **DT**   | ISO8601 date the run was produced                | string |
| **PI**   | Predicted median insert size                     | string |
| **PL**   | Platform or technology used to produce the reads | string |
| **PM**   | Platform model                                   | string |


# Configuration validation

The `--validate` flag makes Pheniqs evaluate the supplied configuration and emit an exhaustive report without actually executing. It is sometimes useful to inspect the report before executing to make sure all implicit parameters are allocated the desired values. The validation report also prints a pairwise Hamming distance matrix for each [multiplex barcode](glossary.html#multiplex barcode) set and one for the concatenated barcode sequences. The top half of the matrix, above the diagonal, is the pairwise Hamming distance, while the bottom half is the maximum number of correctable errors the pair can tolerate. For each barcode set [minimum distance decoder](glossary.html#minimum_distance_decoding) is limited by the smallest value in the bottom half of the matrix.


# Demultiplexing report

Pheniqs emits a comprehensive demultiplexing report with statistics about both inputs and outputs.

## Input statistics

A quality statistics report for every segment in the input is provided in the `demultiplex input report` element.

| JSON field                                      | Description
| :---------------------------------------------- | :----------------------------------------------------------------------
| **count**                                       | Number of input reads
| **pf count**                                    | Number of input reads that [passed vendor quality control](#pass-filter-reads)
| **pf fraction**                                 | **pf count** / **count** 

## Output statistics

A quality statistics report for every segment in every output read group is provided in the `demultiplex output report` element as well as global statistics for the entire pipeline.

### Read Group statistics

Counters in each element of the `read group quality reports` array apply only to reads that were classified to the respective read group.

| JSON field                                      | Description
| :---------------------------------------------- | :---------------------------------------------------------------------
| **count**                                       | Number of reads
| **multiplex distance**                          | Average multiplex distance
| **multiplex confidence**                        | Average multiplex confidence **([pamld](glossary.html#phred_adjusted_maximum_likelihood_decoding) only)**
| **pf count**                                    | Number of reads that [passed vendor quality control](#pass-filter-reads)
| **pf multiplex distance**                       | Average multiplex distance in reads that [passed vendor quality control](#pass-filter-reads)
| **pf multiplex confidence**                     | average multiplex confidence in reads that [passed vendor quality control](#pass-filter-reads) **([pamld](glossary.html#phred_adjusted_maximum_likelihood_decoding) only)**
| **pf fraction**                                 | **pf count** / **count** 
| **pooled fraction**                             | **count** / **pipeline :: count**
| **pf pooled fraction**                          | **pf count** / **pipeline :: pf count**
| **pooled multiplex fraction**                   | **count** / **pipeline :: multiplex count**
| **pf pooled multiplex fraction**                | **pf count** / **pipeline :: pf multiplex count**


### Pipeline statistics

Counters found directly in the `demultiplex output report` element are for the output of the entire pipeline.

| JSON field                                      | Counter incrementing criteria 
| :---------------------------------------------- | :---------------------------------------------------------------------
| **count**                                       | Sum of **count** in all read groups, including undetermined
| **multiplex count**                             | Sum of **count** in all read groups, excluding undetermined
| **multiplex fraction**                          | **multiplex count** / **count**
| **multiplex distance**                          | Average multiplex distance, excluding undetermined
| **multiplex confidence**                        | Average multiplex confidence, excluding undetermined **([pamld](glossary.html#phred_adjusted_maximum_likelihood_decoding) only)**
| **pf count**                                    | Sum of **pf count** in all read groups
| **pf fraction**                                 | **pf count** / **count** 
| **pf multiplex count**                          | Sum of **pf count** in read groups, excluding undetermined
| **pf multiplex fraction**                       | **pf multiplex count** / **pf count**
| **pf multiplex distance**                       | Average multiplex distance in pf reads, excluding undetermined
| **pf multiplex confidence**                     | Average multiplex confidence in pf reads, excluding undetermined **([pamld](glossary.html#phred_adjusted_maximum_likelihood_decoding) only)**
| **multiplex pf fraction**                       | **pf multiplex count** / **multiplex count** 


# Performance tuning

Pheniqs will initialize a thread pool used by [HTSlib](glossary.html#htslib) for IO and an equivalent number or processing threads. Factors affecting the optimal configuration are the length and number of input and output segments, the number of multiplexed libraries, the file format of both input and output files, whether quality tracking has been enabled and the layout of both input and output reads. Reading and writing gzip compressed [FASTQ](glossary.html#fastq) files, for instance, is a bottleneck since gzip is a very old compression algorithm that scales poorly when parallelized. In some scenarios, depending on disk IO speed and space constraints it might be overall more efficient to base call to uncompressed FASTQ files and feed Pheniqs uncompressed FASTQ input. Binary [HTSlib](glossary.html#htslib) formats are generally much faster to encode and decode and should be preferred for both input and output where possible.

| Name                      | Description                                                   | Type    |
| :------------------------ | :------------------------------------------------------------ | :------ |
| `threads`                 | IO thread pool size                                           | integer |
| `transforms`              | number of transforming threads                                | integer |
| `buffer capacity`         | set feed buffer capacity                                      | integer |
| `disable quality control` | do not track quality when processing reads                    | boolean |
