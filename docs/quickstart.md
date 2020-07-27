---
layout: default
title: "Quickstart"
permalink: /quickstart/
id: quickstart
---

* placeholder
{:toc }

The Pheniqs command line interface accepts a [JSON](https://en.wikipedia.org/wiki/JSON) encoded configuration file containing directives for input and output layout, parsing read segments and run parameters. If you are new to JSON, a [validator](cli.html#json-validation) can be instrumental for troubleshooting syntax errors. Some parameters can also be specified as command line arguments that override their corresponding configuration file values. A brief description of the command line parameters Pheniqs accepts is always available with the `-h/--help` flags. If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) the [bundled command line completion](cli.html#zsh-completion) will give you a more interactive command line experience.

# Supported File Format

Pheniqs can arbitrarily manipulate reads from either [SAM, BAM and CRAM](glossary.html#htslib) or [FASTQ](glossary.html#fastq). It can even operate on different read segments in different formats at the same time or interact with segments that are either [interleaved](glossary.html#interleaved_file_layout) into a single file or [split](glossary.html#split_file_layout) over many. Using the intuitive addressing syntax you can classify reads by standard Illumina sample barcodes, multiple combinatorial cellular barcodes and molecular index tags in arbitrary locations along the read, all without pre or post processing the data.

# Declaring Input

In this example we consider three files that contain synchronized segments from an Illumina MiSeq instrument. Pheniqs will construct the input [read](glossary.html#read) by reading one [segment](glossary.html#segment) from each input file. [Relative](glossary.html#relative_path) input and output file paths are resolved against the working directory which defaults to where you execute Pheniqs. You may optionally specify `base input url` to explicitly provide a base directory and make the configuration file more portable.

For example we assume 3 FASTQ files created by executing bcl2fastq to simply get all 3 raw segments of a single indexed Illumina MiSeq run. The first four lines of each of the files are

*000000000-BDGGG_S1_L001_R1_001.fastq.gz*, the first read segment (or *Read 1* in Illumina terminology):
```
@M02455:162:000000000-BDGGG:1:1101:10000:10630 1:N:0:
CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA
+
CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
```

*000000000-BDGGG_S1_L001_I1_001.fastq.gz*, the second read segment (or *Index 1* in Illumina terminology):
```
@M02455:162:000000000-BDGGG:1:1101:10000:10630 2:N:0:
GGACTCCT
+
B@CCCFC<
```

*000000000-BDGGG_S1_L001_R2_001.fastq.gz*, the third read segment (or *Read 2* in Illumina terminology):
```
@M02455:162:000000000-BDGGG:1:1101:10000:10630 3:N:0:
GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA
+
CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
```

To declare those files as input you add an `input` directive, which is a JSON array of file paths. The order of the paths in the array enumerates the input read segments so Pheniqs will expect to read one segment from each file you list. You may provide the input segments in any order you wish as long as you adhere to that order throughout the configuration file.

>```json
{
    "input": [
        "000000000-BDGGG_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_R2_001.fastq.gz"
    ]
}
```
>**Example 1.1** Declaring an input read that is [split](glossary.html#split_file_layout) over three gzip compressed FASTQ files. Since the file paths do not start with `/` they are considered [relative](glossary.html#relative_path) and resolved against the current working directory. If you specify `base input url` paths are resolved relative to that directory path. [Absolute file paths](glossary.html#absolute_path) ignore `base input url`.
{: .example}

Notice that the order of the paths in the array is not just telling Pheniqs where to find the files but actually defines the enumerated segments of the input read. When reading the [split](glossary.html#split_file_layout) read layout in this example that just means Pheniqs will read one segment from each input file. But if the same 3 segment input reads were [interleaved](glossary.html#interleaved_file_layout) into one file, you would list the same file path 3 times to tell Pheniqs that every 3 records in that file form a single sequence read with 3 segments.

**Example 1.1** is already a complete and valid Pheniqs configuration file! Since we have not yet specified any output or manipulation instructions, reads are simply interleaved to the default [stdout](glossary.html#standard_stream) in SAM format. If you execute it you will get the following output on stdout:

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG     FI:i:1  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   13      *       0       0       *       *       0       0       GGACTCCT        B@CCCFC<        FI:i:2  TC:i:3
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG     FI:i:3  TC:i:3
```
>**Example 1.2** Output header and first 3 records (one complete read) from [interleaving](glossary.html#interleaved_file_layout) the three read segments verbatim into a single SAM formatted stream written to standard output using the configuration file in **Example 1.1**.
{: .example}

Next we will show you how to use this simple example to interleave raw split read segments into a single CRAM file. CRAM files are the latest indexed and compressed binary encoding of SAM implemented in [HTSlib](glossary.html#htslib) and often provide more efficient compression than the ubiquitous gzip compressed FASTQ format while being much faster to interact with. Packaging your reads into a CRAM container also makes archiving raw data simple. Another huge advantage of interleaved files is that they may be produced or consumed by Pheniqs through [standard streams](glossary.html#standard_stream).


# Declaring Output

Since most of the time you do not want your output on stdout you will want to provide an output file path. To write interleaved output to a compressed CRAM file simply add an `output` directive to **Example 1.1**. Like `input`, the `output` directive is a JSON array of file paths. You may optionally specify multiple paths in the output array to write a split output.

>```json
{
    "input": [
        "000000000-BDGGG_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_R2_001.fastq.gz"
    ],
    "output": [ "000000000-BDGGG_raw.cram" ]
}
```
>**Example 1.3** Interleaving three read segments verbatim into a single CRAM file. CRAM files are often much faster to read and write, especially in highly parallelized environments, and also support a rich metadata vocabulary.
{: .example}

You may optionally specify `base output url` to explicitly provide a base directory and make the configuration file more portable. [Absolute file paths](glossary.html#absolute_path) ignore `base output url`.

Before we move on to read manipulation here is the same interleaving configuration achieved without a configuration file by specifying the same instruction on the command lines

>```shell
pheniqs mux \
--input BDGGG_s01.fastq.gz \
--input BDGGG_s02.fastq.gz \
--input BDGGG_s03.fastq.gz \
--output 000000000-BDGGG_raw.cram
```
>**Example 1.3.1** Interleaving three read segments verbatim into a single CRAM file without a configuration file. Notice that `--input` is specified 3 times and the order the arguments are provided on the command line enumerates the input segments.
{: .example}

# Read transformation

Pheniqs provides a generic method to derive a new set of sequence segments from the input read. That set can form the desired output or it can be a technical artifact used to classify the biological sequence. Either way, the syntax is the same. When declared inside the `template` section, the `transform` directive constructs the [output](glossary.html#output_segment) read segments. When declared inside a barcode decoder, it constructs the set of sequences that will be assessed against the list of expected barcode sequences.

Transforms operate on the input read in two steps: First the `token` element, a JSON array of [tokenization](manual.html#tokenization) patterns, that each extract a continuous sequence (or a token) from an [input segment](glossary.html#input_segment) and second the `knit` element that [constructs](manual.html#segment-assembly) new segments from the previously defined tokens.

If the segments you are extracting can be found as one continuous sequence in a read segment you only need to specify the `token` array and Pheniqs will assume that each token represents a segment.

If, however, you are trying to extract segments from multiple, non continuous, tokens or need to reverse complement the sequence to match against your expected barcodes you have one more step. The optional [knit](manual.html#transform-pattern) directive references the tokens to construct a new segment from multiple tokens. If the `knit` array is omitted from `transform`, each token is assumed to declare a single segment.

# Output Manipulation

The `transform` directive in the `template` section is used to manipulate the output read. If omitted all segments of the input are written verbatim to the output, as seen in **Example 1.1** and **Example 1.3**. Since the second segment contains only a technical sequence, and we do not want to write it to the output, we add a `transform` directive to construct an output read from only the first and third segments of the input.

>```json
{
    "input": [
        "000000000-BDGGG_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_R2_001.fastq.gz"
    ],
    "template": {
        "transform": { "token": [ "0::", "2::" ] }
    }
}
```
>**Example 1.4** Adding a transform directive composing the output read from only the untouched first and third input segments. Input segments in Pheniqs are indexed and referenced using a [zero based coordinate system](glossary.html#zero_based_coordinate) so the first segment is 0 and the third is 2. Since we want to output the entire first and last segments the start and end coordinates in the pattern are left out to accept their default values.
{: .example}

The [token patterns](manual.html#tokenization) declared in the `token` array of the `transform` directive are made of 3 colon separated integers. The first is the [zero based](glossary.html#zero_based_coordinate) [input segment index](glossary.html#input_segment). The second is an inclusive [zero based](glossary.html#zero_based_coordinate) **start** coordinate to the beginning of the token and it defaults to **0** if omitted. The third is an exclusive [zero based](glossary.html#zero_based_coordinate) **end** coordinate to the end of the token. If the **end** coordinate is omitted the token spans to the end of the segment. The two colons are always mandatory. Pheniqs token pattern can address segments from either the 5' (left) or 3' (right) end. To address the 3' end you use negative coordinates. Since they mimic the [python array slicing](https://en.wikipedia.org/wiki/Array_slicing#1991:_Python) syntax they are fairly easy to test.

>```
@HD     VN:1.0  SO:unknown      GO:query
@RG     ID:undetermined PU:undetermined
M02455:162:000000000-BDGGG:1:1101:10000:10630   77      *       0       0       *       *       0       0       CTAAGAAATAGACCTAGCAGCTAAAAGAGGGTATCCTGAGCCTGTCTCTTA     CCCCCGGGFGGGAFDFGFGGFGFGFGGGGGGGDEFDFFGGFEFGCFEFGEG
M02455:162:000000000-BDGGG:1:1101:10000:10630   141     *       0       0       *       *       0       0       GCTCAGGATACCCTCTTTTAGCTGCTAGGTCTATTTCTTAGCTGTCTCTTA     CCCCCGGGGGGGGGGGGGGGGGGGF<FGGGGGGGGGGGGFGFGGGGGGGGG
```
>**Example 1.5** Output header and first 2 records (one complete read) from interleaving the three read segments to a SAM formatted stream using the configuration file in **Example 1.4**. Notice that only the first and third segments were written to the output.
{: .example}

# Classifying by barcodes

The reads in our files were sequenced from DNA from 5 individually prepared libraries that were tagged with an 8bp technical sequence before they were pooled together for sequencing. To classify them into read groups we need to examine the first 8 nucleotides of the second segment. Decoding the barcode can be as trivial as comparing two strings if we were absolutely confident no errors occurred during sequencing. However in a real world scenario each nucleotide reported by a sequencing instrument is accompanied by an estimate of the probability the base was incorrectly called. We refer to such an uncertain sequence as an **observed sequence**. The `multiplex` directive is used to declare a decoder that will classify the reads by examining the segments constructed by the embedded `transform` and comparing them to the sequence segments declared in the `codec` element.

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
>**Example 1.6** A `multiplex` decoder directive declaration using the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding). The `transform` directive is used to extract observed segments from the raw read while the `codec` directive names the possible barcode sequences we expect to find.
{: .example}

In this example we declare a `multiplex` directive that uses the [phred-adjusted maximum likelihood decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding) algorithm. This algorithm will choose a barcode using a maximum likelihood estimate and reject any classification with a decoding confidence lower than the `confidence threshold` parameter. The `noise` parameter is the prior probability that an observed sequence has not originated from any of the provided barcodes. The value of `noise` is often set to the amount of [PhiX Control Library](http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf) spiked into the solution for reads sequenced on the Illumina platform but can be higher if you expect other types of noise to be present. In the `codec` directive we provide a discrete set of possible decoding results. All `barcode` segment arrays must match the layout declared in the embedded `transform` directive, in this example one 8bp segment. The keys of the `codec` directive can be any unique string. In this example we used the unique barcode nucleotide sequence prefixed with an @ character to remind us this is simply a unique identifier.

>```json
{
    "input": [
        "000000000-BDGGG_S1_L001_R1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_I1_001.fastq.gz",
        "000000000-BDGGG_S1_L001_R2_001.fastq.gz"
    ],
    "template": {
        "transform": { "token": [ "0::", "2::" ] }    
    },
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
    "template": {
        "transform": { "token": [ "0::", "2::" ] }
    },
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
