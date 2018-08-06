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
        <li><a                  href="/pheniqs/2.0/best_practices.html">Tips and Best Practices</a></li>
        <li><a                  href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Build</a></li>
        <li><a                  href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a class="active"   href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

<section class="notice_2_0">Please excuse us as we update the documentation for the new Pheniqs 2.0 API.</section>

# Pheniqs Documentation
{:.page-title}

* placeholder
{:toc}

# Configuration
For most non trivial scenarios Pheniqs requires a [JSON](https://en.wikipedia.org/wiki/JSON) encoded instruction file. Some parameters you can interactively override with [command line](cli.html) arguments and almost all instruction directives have a default value. The instruction file contains a number of separate sections specifying directives for declaring input and output layout, read segment manipulation, barcode decoding and run parameters. Since barcode decoding often relies on verbose configuration semantics, Pheniqs supports a sophisticated instruction inheritance mechanism and can also import additional instruction files when compiling your job. While you may certainly ignore this added complexity at first, separating different aspects of your experimental design into reusable instruction files can significantly simplify day to day use.

# Format support and performance
Pheniqs can arbitrarily manipulate reads stored in either the advanced [HTS](glossary.html#htslib) containers SAM, BAM and CRAM or the legacy [FASTQ](glossary.html#fastq) format with segments either [interleaved](glossary.html#interleaved_file_layout) into a single file or [split](glossary.html#split_file_layout) over many. Some features are only supported with HTS containers, since FASTQ records have no standardized method of associating metadata with read segments. Pheniqs is written in highly optimized multi threaded [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) and interfaces directly with the low level [HTSlib](glossary.html#htslib) API. It can be installed from many popular package managers or compiled from source into a single portable executable binary using the provided `pheniqs-tools` python script.

Pheniqs execution speed will normally scale linearly with the number of cores available for computation. However since the ancient gzip compression algorithm does not scale well in multi threaded environments, reading or writing from gzip compressed FASTQ files will often be I/O bound. You can use Pheniqs to easily repack legacy data stored in FASTQ files into an efficient, interleaved, CRAM container. CRAM offers significant improvements in both performance and storage requirements as well as support for associating extensible metadata with read segments.

Pheniqs aims to fill the gap between advanced HTS containers and legacy analysis tools that only supports FASTQ. By efficiently converting between HTS and FASTQ over [standard streams](https://en.wikipedia.org/wiki/Standard_streams), it allows you to feed legacy analysis software with FASTQ directly from annotated HTS files, without additional storage requirements.

# The `input` directive
The instruction `input` directive is an ordered list of file paths. Pheniqs assembles an input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input file.

>```json
{
    "input": [
        "HK5NHBGXX_Lane1_S1_L001_R1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I2_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_R2_001.fastq.gz",
    ]
}
```
>**Example 2.1** Declaring an input read that is [split](glossary.html#split_file_layout) over four gzip compressed FASTQ files.
{: .example}

[Interleaved](glossary.html#interleaved_file_layout) files contain multiple consecutive segments of the same read. To explicitly assemble an input read from interleaved files simply repeat the path to reference the same file multiple times, once for each [segment](glossary.html#segment).

>```json
{
    "input": [
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram"
    ]
}
```
>**Example 2.2** Declaring a 4 segment input read [interleaved](glossary.html#split_file_layout) into a single CRAM file.
{: .example}

You can even mix-and-match the two layout styles

>```json
{
    "input": [
        "HK5NHBGXX_Lane1_biological.fastq.gz",
        "HK5NHBGXX_Lane1_technical.fastq.gz",
        "HK5NHBGXX_Lane1_technical.fastq.gz",
        "HK5NHBGXX_Lane1_biological.fastq.gz",
    ]
}
```
>**Example 2.3** Constructing a four segment  input read from a mixed style layout. The two biological segments, 0 and 3, are [interleaved](glossary.html#interleaved_file_layout) in `HK5NHBGXX_Lane1_biological.fastq.gz` while the two technical segments, 1 and 2, are interleaved into `HK5NHBGXX_Lane1_technical.fastq.gz`.
{: .example}

The `input` directive defaults to expecting an interleaved SAM file provided on standard input

>```json
{
    "input": [ "/dev/stdin" ]
}
```
>**Example 2.4** Pheniqs implicitly expects input on standard input if not instructed otherwise. The format and layout will be automatically detected.
{: .example}

## Automatic `input` sensing
Pheniqs will automatically detect the format of the input you provide it by examining the first few bytes. Once the input format is established, Pheniqs decodes the first few segment records from the feed to detect the feed's resolution, or the number of consecutive segments in the feed that have the same read identifier. Pheniqs will assume that to assemble a read it will have to read that many segments from each input feed.

>```
@M02455:162:000000000-BDGGG:1:1101:10000:10630 1:N:0:
CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA
+
CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
@M02455:162:000000000-BDGGG:1:1101:10000:10630 2:N:0:
GGACTCCT
+
B@CCCFC<
@M02455:162:000000000-BDGGG:1:1101:10000:10630 3:N:0:
GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA
+
CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
@M02455:162:000000000-BDGGG:1:1101:10000:12232 1:N:0:
GTATAGGGGTCACATATAGTTGGTGTGCTTTGTGAACTGCGATCTTGACGG
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@M02455:162:000000000-BDGGG:1:1101:10000:12232 2:N:0:
GGACTCCT
+
CCCCCGGG
@M02455:162:000000000-BDGGG:1:1101:10000:12232 3:N:0:
GTCCTATCCTACTCGGCTTCTCCCCATTTTTCAGACATTTTCCTATCAGTC
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
```
>**Example 2.5** This interleaved FASTQ file has a resolution value of **3**.
{: .example}

# The `transform` directive
Read manipulation is achieved with the `transform` directive by means of [tokenization](#tokenization) with the embedded `token` directive and [segment assembly](#segment_assembly) with the `segment pattern` directive. In the tokenization step Pheniqs consults a token pattern to extract a token from an [input segment](glossary.html#input_segment). In the segment assembly step one or more [segment patterns](manual.html#transform-pattern) reference tokens to assemble a new segment. The optional `segment pattern` directive is only necessary when assembling segments from multiple, non contiguous, tokens or if a token needs to be reverse complemented. If the `segment pattern` directive is omitted from a `transform` directive, each token implicitly declares a single segment.

## Tokenization
A token pattern is made of 3 colon separated integers. The first is the mandatory [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment) enumerated by `input`. The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. **start** coordinate and **end** coordinate can take positive or negative values to access the segment from either the 5' (left) or 3' (right) end and mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax. The two colons are always mandatory.

>```json
{
    "transform": { "token": [ "0::", "1::8" ] }
}
```
>**Example 2.6** In this `transform` directive the first token spans the entire first input segment. The second spans the first 8 cycles of the second segment. Since the `segment pattern` directive is omitted each token effectively declares a segment.
{: .example}

>```json
{
    "transform": {
        "token": [ "0:0:", "1:0:8" ],
        "segment pattern": [ "0", "1" ]
    }
}
```
>**Example 2.7** This `transform` directive is semantically identical to the one provided in **Example 2.6** but explicitly declares the `segment pattern` directive and token **start** coordinate.
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
>**Example 2.7** Some examples of tokenizing a hypothetical segment with index 0 and 9 cycles.
{: .example}

## Segment assembly
A `segment pattern` is made of one or more token references separated by the **:** concatenation operator. A token reference is the [zero based](glossary.html#zero_based_coordinate) index of the token pattern in the adjacent `token` directive. Appending the left hand side reverse complementarity **~** operator will concatenate the reverse complemented sequence of the token. Each token reference is evaluated before concatenation so **~** evaluation precedes **:** evaluation.

>| Pattern | Description                                                                                        |
>| ------- | :------------------------------------------------------------------------------------------------- |
>| `0`     | The simplest possible transform will assemble a segment from just token 0.                         |
>| `~0`    | Assemble a segment from the reverse complemented token 0.                                          |
>| `0:1`   | Assemble a segment by concatenating token 0 and token 1                                            |
>| `~0:1`  | Assemble a segment by concatenating the reverse complement of token 0 and token 1.                 |
>
>**Example 2.8** Several examples of `segment pattern` syntax for constructing new segments.
{: .example}

Notice that a negative token **start** coordinate is equivalent to the corresponding positive token **end** coordinate on the reverse complemented strand and vice verse. More formally the token `0:-x:-y` is equivalent to `0:y:x` if applied to the reverse complement of the segment. For instance to concatenate to the first output segment the first 6 bases of the reverse complemented strand of the first input segment you would define token `0:-6:` and then reference it in the transform pattern for output segment 0 as `~0`.

## Contextual `transform` directives
A `transform` directive declared in the root of the instruction assembles the [output read](glossary.html#output_segment) and implicitly defines the `output segment cardinality`. When declared inside a [barcode decoder](#barcode-decoding) directive a `transform` constructs the segmented sequence that is matched against the possible barcode values.

>```json
{
    "input": [
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram"
    ],
    "transform": {
        "token": [ "0:6:", "3::-6" ]
    },
    "multiplex": {
        "transform": { "token": [ "1::8", "2::8" ] }
    },
    "molecular": [
        {
            "transform": {
                "token": [ "0::6", "3:-6:" ],
                "segment pattern": [ "4", "~5" ]
            }
        }
    ]
}
```
>**Example 2.9** Consider that the input we declared in **Example 2.2**, is a dual indexed paired end read with two 8 cycle multiplex barcode segments on the second and third input segments and two molecular barcode segments, one on the first 6 cycles of the first segment and another, reverse complemented, on the last 6 cycles of the last segment. We extract both molecular barcode segments and remove them from the output segments.
{: .example}

## The `output` directive
The instruction `output` directive is an ordered list of output file paths. If the output directive contains the same number of file paths as the `output segment cardinality` each segment will be written to the corresponding enumerated output file. If only one file path is declared all segments will be interleaved into that file. Declaring a number of files paths different than **1** or `output segment cardinality` is ambiguous and will result in a validation error. Read segments are guaranteed to be written contiguously and in-order. The order of the reads in the output is **not** guaranteed to be consistent with their order in the input or to be stable between successive executions when the `thread` directive is bigger than **1**. The `output` directive defaults to interleaved SAM written to standard output.

# Barcode decoding
Pheniqs can populate the related standardized SAM auxiliary tags for multiplex, molecular and cellular barcodes, adhering to the [recommendations outlined in the SAM specification](https://samtools.github.io/hts-specs/SAMtags.pdf). We distinguish between [closed class decoding algorithms](glossary.html#closed_class_decoding), when a discrete list of classes (in our case nucleotide sequences) is known in advance and so a prior distribution is available, and [open class decoders](glossary.html#open_class_decoding), when the discrete list of classes is unknown and so the prior distribution is hidden. The embedded `codec` directive specifies the allowed classes for closed class decoders as a JSON dictionary. The keys of the dictionary must be unique strings within the dictionary and play a role in the inheritance model but you may otherwise choose them as you see fit. One simple methodology is to use the concatenated barcode sequence prefixed by an **@** sign (to remind you Pheniqs does not actually interpret it as a nucleotide sequence), but you may choose more meaningful names to make your instruction files more readable.

Pheniqs offers a choice of two closed class decoding strategies: the widespread [minimum distance decoder](glossary.html#minimum_distance_decoding) (**MDD**) and a more refined probabilistic [Phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) (**PAMLD**). Unlike MDD, PAMLD consults base calling quality scores and the prior barcode distribution and will outperform MDD in almost every real world scenario. Since PAMLD computes the full Bayesian decoding probability to pick the maximum likelihood, the probability of an incorrect barcode assignment is made available in SAM auxiliary tags reserved for local use, for downstream analysis consideration.

Pheniqs currently does not support degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) bases in closed class barcode declarations, but the probabilistic model underlaying PAMLD can be extended to support them if a demand arises.

Since error handling is much trickier in open class decoding, Pheniqs currently only supports a [naive](glossary.html#naive_decoding) decoder that will simply populate the relevant raw, uncorrected, SAM auxiliary tags without attempting error correction. In the case of decoding open class molecular barcodes, all reads associated with a unique molecular identifier are expected to be tagging the same DNA molecule. For that reason, most published error correction strategies involve examining the corresponding biological sequence and aligning it to a consensus sequence of a group of reads that have already been associated with that barcode. A specialized open class error correcting molecular barcode decoder is planned for a future release.

When declaring a decoder directive you specify the decoding strategy by setting the `algorithm` attribute.

## Handling decoding failure
The `undetermined` directive can be used to inform a decoder how to handle reads that fail decoding. If omitted the implicit `undetermined` directive inherits global or decoder attributes in the same way that specific barcode classes do. When populating corrected nucleotide and quality sequences in SAM auxiliary tags, a placeholder sequence of **=** character of the correct layout is used to denote a decoding failure.

## Phred-adjusted maximum likelihood decoding
Pheniqs implements a probabilistic closed class [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) that directly estimates the decoding likelihood from the base calling error probabilities provided by the sequencing platform. When output is SAM encoded Pheniqs will write the decoding error probability to the [XB](glossary.html#dq_auxiliary_tag) auxiliary tag for a multiplex barcode,  the [XM](glossary.html#xm_auxiliary_tag) tag for a molecular barcode and the [XC](glossary.html#xc_auxiliary_tag) tag for a cellular barcode. In the case of a decoding failure no error probability is reported.

The dominant parameter to set when decoding with PAMLD is the `confidence threshold`, which is a lower bound on the decoding probability for the decoder to declare a successful classification. The value of `confidence threshold` defaults to **0.99** but in practice depends on the application. Since the decoding error probability is written to the a SAM auxiliary tag you may choose to set `confidence threshold` to a relatively lax value and leave the decision to exclude certain reads for downstream analysis.

PAMLD also consults user provided prior barcode instance and decoding failure frequencies. Barcode frequencies are specified in the class `concentration` attribute while the expected decoding failure frequency is specified in the decoder `noise` attribute. For the Illumina platform the most common value for `noise` will be the amount of [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the pooled solution, which is usually **1%** and corresponds to a value of **0.01**. Assays with low base diversity sometimes use a higher concentration to compensate for low base diversity and some applications can go as high as **50%**. If you know a priori that some amount of sequences present in the pool represent contamination or multiplexed libraries you are not declaring in the configuration, their concentration should be factored into this parameter.

The `concentration` attribute defaults to **1** if omitted. The values for all barcode instances in a decoder are normalized so that their sum equals **1.0** - `noise`. Notice that unlike `concentration` the `noise` attribute is specified as a probability value between **0** and **1**, and it rarely make sense to set it to **0**. If the `concentration` attribute is omitted from all classes the result is an implicit uniformly distributed prior.

## Minimum distance decoding
Pheniqs also implements the more traditional discrete [minimum distance decoder](glossary.html#minimum_distance_decoding) that consults the edit distance between the expected and observed sequence. MDD consults the `distance tolerance` attribute, which is a list of upper bounds on the edit distance between each segment of the expected and observed barcode to still be considered a match. Setting this property to a value higher than the [Shannon bound](glossary.html#shannon_bound), which also serves as the default value for `distance tolerance`, can lead to ambiguous classification and will result in a validation error.

Since MDD effectively ignores the Phred encoded quality scores, it may be consulting extremely unreliable base calls. To mitigate that effect you may set the `quality masking threshold` attribute, which is a lower bound on the permissible base calling quality. Observed bases with quality lower than this threshold will be considered as **N** by the minimum distance decoder. `quality masking threshold` defaults to **0** which effectively disables quality masking.

## The `multiplex` directive
A single closed class decoder can be declared in the `multiplex` directive. When decoding multiplex barcodes Pheniqs will write the nucleotide barcode sequence to the [BC](glossary.html#bc_auxiliary_tag) SAM auxiliary tag, the corresponding Phred encoded quality sequence to the [QT](glossary.html#qt_auxiliary_tag) tag and, when decoding with PAMLD, the decoding error probability to the [XB](glossary.html#xb_auxiliary_tag) tag. Pheniqs classifies the read to a read group by populating the [RG](glossary.html#rg_auxiliary_tag) SAM auxiliary tag, which is a reference to the **ID** attribute of a read group declared in the SAM header.

>```json
{
    "input": [
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram"
    ],
    "multiplex": {
        "transform": { "token": [ "1::8", "2::8" ] },
        "codec": {
            "@AAGAGGCAAGAGGATA": { "barcode": [ "AAGAGGCA", "AGAGGATA" ] },
            "@AGGCAGAAAGAGGATA": { "barcode": [ "AGGCAGAA", "AGAGGATA" ] },
            "@CAGAGAGGAGAGGATA": { "barcode": [ "CAGAGAGG", "AGAGGATA" ] },
            "@CGTACTAGTCTTACGC": { "barcode": [ "CGTACTAG", "TCTTACGC" ] }
        }
    }
}
```
>**Example 2.10** Expanding **Example 2.9** we declare 4 classes for the multiplex decoder to be matched against the segmented sequence extracted by the `transform`. The keys of the `codec` dictionary directive have no special meaning and you may choose them as you see fit, as long as they are unique within the `codec` dictionary.
{: .example}

Since each class decoded by the multiplex decoder corresponds to a read group you may define read group related attributes in each class entry in the `codec` dictionary. Read group attributes that apply to all read groups may be declared upstream, in either the multiplex decoder directive or the root of the instruction document.

>```json
{
    "flowcell id": "HK5NHBGXX",
    "flowcell lane number": 1,
    "input": [
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram",
        "HK5NHBGXX_Lane1.cram"
    ],
    "multiplex": {
        "CN": "NYU CGSB",
        "DT": "2016-07-15T07:00:00+00:00",
        "PI": "500",
        "PL": "ILLUMINA",
        "PM": "HiSeq 2500",
        "transform": { "token": [ "1::8", "2::8" ] },
        "codec": {
            "@AAGAGGCAAGAGGATA": {
                "barcode": [ "AAGAGGCA", "AGAGGATA" ],
                "SM": "c57bl6 mouse 1",
                "LB": "pre b fetal liver technical replicant 1"
            },
            "@AGGCAGAAAGAGGATA": {
                "barcode": [ "AGGCAGAA", "AGAGGATA" ],
                "SM": "c57bl6 mouse 1",
                "LB": "pre b fetal liver technical replicant 2"
            },
            "@CAGAGAGGAGAGGATA": {
                "barcode": [ "CAGAGAGG", "AGAGGATA" ],
                "SM": "c57bl6 mouse 1",
                "LB": "pre b fetal liver technical replicant 3"
            },
            "@CGTACTAGTCTTACGC": {
                "barcode": [ "CGTACTAG", "TCTTACGC" ],
                "SM": "c57bl6 mouse 2",
                "LB": "follicular spleen technical replicant 1",
                "PI": "300"
            }
        }
    }
}
```
>**Example 2.11** Further expanding **Example 2.10** with attributes related the read group SAM header tags. Attributes declared in the decoder directive will apply to all barcode entries in the `codec` dictionary unless explicitly overridden in the entry. For instance all except **@CGTACTAGTCTTACGC**, that locally overrides **PI** to be **300**, will have their **PI** tag set to **500**.
{: .example}

| Name                                      | Description                                      | Type   |
| :---------------------------------------- | :----------------------------------------------- | :----- |
| **ID**                                    | Read group identifier                            | string |
| **[LB](glossary.html#lb_auxiliary_tag)**  | Library name                                     | string |
| **[SM](glossary.html#sm_auxiliary_tag)**  | Sample name                                      | string |
| **[PU](glossary.html#pu_auxiliary_tag)**  | Platform unit unique identifier                  | string |
| **CN**                                    | Name of sequencing center producing the read     | string |
| **DS**                                    | Description                                      | string |
| **DT**                                    | [ISO8601](https://en.wikipedia.org/wiki/ISO_8601) date the run was produced                | string |
| **PI**                                    | Predicted median insert size                     | string |
| **[PL](glossary.html#pl_auxiliary_tag)**  | Platform or technology used to produce the reads | string |
| **PM**                                    | Platform model                                   | string |
| **PG**                                    | Programs used for processing the read group      | string |

* **ID** defaults to the value of **PU** if not explicitly specified. If explicitly declared it must be unique within the `codec` dictionary.
* **PU**, following the [convention established by GATK](https://software.broadinstitute.org/gatk/guide/article?id=6472), defaults to `flowcell id`:`flowcell lane number`: `barcode`. If `flowcell id` or `flowcell lane number` are not specified they are omitted along with their trailing `:`. If explicitly declared **PU** must be unique within the `codec` dictionary.
* **PL**, as defined in the SAM specification, is one of *CAPILLARY*, *LS454*, *ILLUMINA*, *SOLID*, *HELICOS*, *IONTORRENT*, *ONT*, *PACBIO*.


## Contextual `output` directives
The `output` attribute can be declared in the root of the instruction, in the `multiplex` decoder directive or in each of the individual `multiplex` decoder class directives. As always, attributes declared deeper in the hierarchy override attributes declared upstream. When interleaving reads from multiple read groups into the same output it is sufficient, and less verbose, to declare the output attribute upstream. When splitting reads from different read groups to different output files you declare the `output` attribute individually for the read group. You may also mix-and-match the two styles. A corresponding [@RG header tag](glossary.html#rg_header_tag) will be added to the header of an output file if at least one segment of reads tagged with that read group are written to it.

## The `molecular` directive
The `molecular` directive can be used to declare either a single decoder or and array containing multiple decoders. When decoding molecular barcodes Pheniqs will write the nucleotide barcode sequence to the [RX](glossary.html#rx_auxiliary_tag) SAM auxiliary tag and the corresponding Phred encoded quality sequence to the [OX](glossary.html#ox_auxiliary_tag) tag. A molecular identifier is written to the [MI](glossary.html#mi_auxiliary_tag) tag.

If error correction is applied to the barcode the raw, uncorrected, nucleotide and quality sequences are written to the [QX](glossary.html#qx_auxiliary_tag) and [BZ](glossary.html#bz_auxiliary_tag) tags and the decoding error probability to the [XM](glossary.html#xm_auxiliary_tag) tag.

## The `cellular` directive
The `cellular` directive can be used to declare either a single decoder or an array containing multiple decoders. When decoding cellular barcodes Pheniqs will write the raw, uncorrected, nucleotide barcode sequence to the [CR](glossary.html#cr_auxiliary_tag) SAM auxiliary tag and the corresponding Phred encoded quality sequence to the [CY](glossary.html#cy_auxiliary_tag) tag, while The decoded cellular barcode is written to the [CB](glossary.html#cb_auxiliary_tag) tag. The decoding error probability is written to the [XC](glossary.html#cr_auxiliary_tag) tag.

# URL handling
Setting global URL prefixes make your instruction file more portable. If specified, the `base input url` and `base output url` are used as a prefix to **relative** URLs defined in the `input` and `output` directives respectively. A URL is considered relative if it **does not** begin with a **/** character. Environment variables in URLs will be resolved by Pheniqs when it compiles your instruction file. `base input url` and `base output url` default to the `working directory` which is the directory where pheniqs was executed. **relative** URLs are resolved against the `working directory`.

>```json
{
    "base input path": "/volume/alpha/HK5NHBGXX/lane_01",
    "input": [
        "HK5NHBGXX_Lane1_S1_L001_R1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I2_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_R2_001.fastq.gz",
    ]
}
```
>**Example 2.12** Providing a `base input path` and `base output path` and **relative** URLs in the `input` or `output` directives is a good way to make your instruction file more portable.
{: .example}

>```json
{
    "input": [
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_Lane1_S1_L001_R1_001.fastq.gz",
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_Lane1_S1_L001_I1_001.fastq.gz",
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_Lane1_S1_L001_I2_001.fastq.gz",
        "/volume/alpha/HK5NHBGXX/lane_01/HK5NHBGXX_Lane1_S1_L001_R2_001.fastq.gz",
    ]
}
```
>**Example 2.13** Compiling the directives in **Example 2.12**.
{: .example}

>```json
{
    "base input path": "~/HK5NHBGXX",
    "input": [
        "HK5NHBGXX_Lane1_S1_L001_R1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I1_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_I2_001.fastq.gz",
        "HK5NHBGXX_Lane1_S1_L001_R2_001.fastq.gz",
    ]
}
```
>**Example 2.14** `base input path` and `base output path` do not have to be **absolute** URLs and may contain environment variables enclosed in curly brackets. The `~` character in the beginning of a URL is interpreted as the `HOME` environment variable on most POSIX shells and resolves to the home directory of the currently logged in user. **relative** URLs are resolved against the `working directory` which is the directory where pheniqs was executed.
{: .example}

## Standard streams
If the `input` directive is omitted Pheniqs will expect input on **/dev/stdin**. If the `output` directive is omitted Pheniqs will write to **/dev/stdout**. Regardless of file extensions, Pheniqs will detect input format by inspecting the first few bytes of your input and will attempt to guess the input resolution and layout. Since Pheniqs selects the output format by inspecting the extensions of the URLs you provide in the `output` directive, specifying output as **/dev/stdout** or **/dev/stderr** does not allow you to specify an output format. In this case you can use the **exploded** URL directive syntax to provide more details about the output stream.

>```json
{
    "output": [
        { "path": "/dev/stdout", "type": "cram" }
    ]
}
```
>**Example 2.15** Declaring interleaved CRAM output to standard output using the **exploded** URL syntax.
{: .example}

# Phred offset
The `input phred offset` and `output phred offset` are applicable only to [FASTQ](glossary.html#fastq) files and specify the Phred scale decoding and encoding [offset](https://en.wikipedia.org/wiki/FASTQ_format#Encoding), respectively. The default value for both is **33**, knowns as the [Sanger format](glossary.html#sanger_format). The binary [HTSlib](glossary.html#htslib) formats BAM and CRAM encode the quality value numerically and so require no further manipulation. The [sequence alignment map format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) states that text encoded SAM records always use the Sanger format.

>```json
{
    "output phred offset": 33,
    "input phred offset": 33
}
```
> Both the `input phred offset` and `output phred offset` default to 33, known as the Sanger format.
{: .example}

# Leading Segment
Pheniqs constructs new segments for the output read and must generate corresponding identifiers and replicate metadata from the input to the output segments. Since segments can potentially disagree on metadata, one input segment is elected as the leader and used as a template when constructing the output segments. The leading segment property is a reference to an [input segment index](glossary.html#input_segment) and defaults to **0** if omitted.

>```json
{ "leading segment": 0 }
```
>The leading segment defaults to be the first input segment.
{: .example}

# Pass filter reads
Some reads are marked by the sequencing platform as not passing the vendor quality control. For instance Illumina sequencers perform an internal [quality filtering procedure](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-percent-pf-technical-note-770-2014-043.pdf) called chastity filter, and reads that pass this filter are called **PF** for **pass-filter**. This can be [signaled](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) on the comment portion of the read identifier in [FASTQ](glossary.html#fastq) files or the **512** flag on a SAM record flag. Setting this attribute to **true** will instruct Pheniqs to decode those unfiltered reads and include them in the output. If the attribute is set to **false** those reads are dropped. `include filtered` defaults to **false**.

# Configuration validation
The `-V/--validate` command line flag makes Pheniqs evaluate the supplied instruction and emit a human readable description of the instruction without actually executing it. It is sometimes useful to inspect this description before executing to make sure all implicit parameters are allocated the desired values. To also print out the barcode distance metric for each closed class decoder you may additionally set the `-D/--distance` command line flag. The top half of the matrix, above the diagonal, is the pairwise Hamming distance, while the bottom half is the maximum number of correctable errors the pair can tolerate, known as the Shannon bound.

# Demultiplexing statistics
Pheniqs emits a comprehensive demultiplexing report with statistics about both inputs and outputs.

## Input
A quality statistics report for every segment in the input is provided in the `demultiplex input report` element.

| JSON field                                      | Description
| : --------------------------------------------- | :---------------------------------------------------------------------
| **count**                                       | all input reads
| **pf count**                                    | input reads that *passed vendor quality control*
| **pf fraction**                                 | **pf count** / **count**

## Output
A quality statistics report for every segment in every output read group is provided in the `demultiplex output report` element as well as global statistics for the entire pipeline.

### Read Group
Counters in each element of the `read group quality reports` array apply only to reads that were classified to the respective read group.

| JSON field                                      | Description
| : --------------------------------------------- | :---------------------------------------------------------------------
| **count**                                       | all reads classified to the read group
| **multiplex distance**                          | average multiplex distance
| **multiplex confidence**                        | average multiplex confidence
| **pf count**                                    | reads that *passed vendor quality control*
| **pf multiplex distance**                       | average multiplex distance in reads that *passed vendor quality control*
| **pf multiplex confidence**                     | average multiplex confidence in reads that *passed vendor quality control*
| **pf fraction**                                 | **pf count** / **count**
| **pooled fraction**                             | **count** / **pipeline :: count**
| **pf pooled fraction**                          | **pf count** / **pipeline :: pf count**
| **pooled multiplex fraction**                   | **count** / **pipeline :: multiplex count**
| **pf pooled multiplex fraction**                | **pf count** / **pipeline :: pf multiplex count**

### Pipeline statistics
Counters found directly in the `demultiplex output report` element are for the output of the entire pipeline.

| JSON field                                      | Counter incrementing criteria
| :---------------------------------------------- | :---------------------------------------------------------------------
| **count**                                       | sum of **count** in all read groups
| **multiplex count**                             | sum of **count** in all read groups, excluding undetermined.
| **multiplex fraction**                          | **multiplex count** / **count**
| **multiplex distance**                          | average multiplex distance, excluding undetermined.
| **multiplex confidence**                        | average multiplex confidence, excluding undetermined.
| **pf count**                                    | sum of **pf count** in all read groups.
| **pf fraction**                                 | **pf count** / **count**
| **pf multiplex count**                          | sum of **pf count** in read groups, excluding undetermined.
| **pf multiplex fraction**                       | **pf multiplex count** / **pf count**
| **pf multiplex distance**                       | average multiplex distance in pf reads, excluding undetermined.
| **pf multiplex confidence**                     | average multiplex confidence in pf reads, excluding undetermined.
| **multiplex pf fraction**                       | **pf multiplex count** / **multiplex count**
