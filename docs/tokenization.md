---
layout: default
title: "Tokens and Transforms"
permalink: /tokenization/
id: tokenization
---

# Barcode Tokens and Transform Patterns
{:.page-title}

* placeholder
{:toc}

Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experimental designs.

The conceptual framework that Pheniqs uses to extract different types of elements within sequence reads is summarized below.

# Experimental Design

Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. Some very common designs for Illumina platforms are illustrated here:

![experimental designs](/pheniqs/assets/img/diagram8.png)

More complicated barcoding schemes for single-cell, CRISPR, and multi-modal sequencing are also appearing. Examples of how to configure Pheniqs for a variety of experimental designs may be found in the [vignettes section](vignettes.html) of the documentation.


# Read anatomy

Illumina sequencing platforms typically produce four different sequence elements: two Index sequences, referred by Illumina as the **i5** and **i7** barcodes, and two Insert sequences, referred by Illumina as **read 1** and **read 2**. Collectively, these are referred to as read segments. For example, consider a [standard paired-end, dual index library design](illumina.html):

![read anatomy](/pheniqs/assets/img/diagram1.png)

The read segments for this standard design thus comprise two technical sequences (referred by Illumina as I1, I2) and two biological sequences (referred by Illumina as R1, R2):

![read anatomy](/pheniqs/assets/img/diagram2.png)


## Sequence Classification

The combination of the barcodes contained in the **I1** and **I2** index positions specifies the sample library. With standard dual indexing, up to 96 distinct sample libraries can be pooled and run together in a single sequencing lane.

To identify which biological sequences belong to which library, the sequences belonging each one need to be separated from each other. This process of deconvolving libraries is called demultiplexing and is done by classifying each of the sequences using the barcode indexes:

![read anatomy](/pheniqs/assets/img/diagram5.png)
>**Note** The i5 adaptor sequences specified in sample sheets will be reverse complemented for platforms that read the I2 index on the bottom strand. The i7 sequences in sample sheets are always reverse complemented relative to the original adaptor sequences since they are read from the top strand.
{: .example}


## Tokenization

Pheniqs uses [tokens](manual.html#tokenization) to reference and extract information from different read segments by specifying where to look for different classes of sequence elements (i.e. barcodes, biological sequences). Each element of interest is defined by an offset relative to the beginning of a given read segment (in this example I1, I2, R1, R2) and a length. It is important to note that Pheniqs uses [zero based](glossary.html#zero_based_coordinate) indexing, so the first read to come off the machine will be referred to as Segment 0, and so on:

![read anatomy](/pheniqs/assets/img/diagram7.png)

>For this design, the barcode tokens begin at position 0 in I1 and I2 and extend for 10 bases. The tokens for biological sequences begin at position 0 of Read 1 and Read 2 and extend for the full span of those read segments.
{: .example}

For a standard paired-end, dual indexed Illumina run, the sample barcodes usually comprise the full I1 and I2 read segments. Because Illumina sequencing is calibrated in relation to the previously sequenced base, those segments are sometimes sequenced one nucleotide longer than necessary to ensure good quality on the last nucleotide. The biological sequences start at the first position of R1 and R2 and extend for the full number of cycles run (typically 75, 100, or 150 nucleotides).


## Transform Patterns

This example illustrates tokenization syntax and output for a 150nt dual-indexed paired-end sequencing run with sample, cellular, and molecular barcodes. This example contains the following features:

+ A **sample** barcode composed of two 10nt elements (i5 and i7)
+ A 12nt inline **cellular** barcode (Cell)
+ A 12nt inline **molecular** barcode (UMI)
+ An **Insert** containing the biological sequence of interest (**template**), which is sequenced from both ends. The template sequences are located in Read 1 (31nt just downstream of the Cell and UMI) and all of Read 2 (here, 75nt).

![transform patterns](/pheniqs/assets/img/transforms.png)

>_**Input files**_ contain sequence data to be analyzed. For Illumina data, these will usually be in fastq format and there will be one file for each read _**segment**_ that comes off the sequencing machine.
The read _**segments**_ are catalogued for future reference and are indexed as an array, where 0=Read1, 1=Index1, 2=Index2, 3=Read2.
Barcode _**tokens**_ define the type and location of the sequences of interest in each read segment. Barcodes may appear at any position and orientation in any read segment. Each type of barcode included in the experimental design is designated as a separate entity.
The _**transform patterns**_ define how each sequence component to be extracted is to be handled. Each token comprises three colon separated components, ``segment:start:end``. Per Python array slicing syntax, the *start* coordinate (offset) is inclusive and the *end* coordinate is exclusive. Start and end coordinates default to 0 and the end of the segment, respectively.
{: .example}

## Output

After sequence classification, the tokenized sequences that have been extracted are written along with their associated metadata, including sequence quality and confidence scores, to specific [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf) within the [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format. (While Pheniqs can also produce FASTQ output, SAM is preferred since it preserves all associated metadata for each read group.)

<img src="/pheniqs/assets/img/sam_output.png" width="400" />

>Sample output for the 150nt paired-end dual-indexed experimental design shown above. Template read segments are emitted along with observed and most likely inferred barcode sequences, quality scores, and error probabilities. The _**confidence score**_ for each token is one minus its estimated error based on the full posterior probability of observation; for the compound sample barcode here, it is the product of the confidence scores for each component and is one minus the error probability shown.
{: .example}

The SAM tags populated by Pheniqs are summarized below:

| Name                                      | Description                                                    | Example                       |
| :---------------------------------------- | :------------------------------------------------------------- | :---------------------------- |
| **[RG](glossary.html#rg_auxiliary_tag)**  | Read group identifier matching an RG filed in the header.      | H7LT2DSXX:1:GAACTGAGCGTCGTGGAGCG  |
| **[BC](glossary.html#bc_auxiliary_tag)**  | Raw uncorrected sample barcode sequence.                       | GAACTGAGCG-TCGTGGAGCG             |
| **[QT](glossary.html#qt_auxiliary_tag)**  | Phred quality of the sample barcode sequence in the BC tag.    | ,FF::F:F:F-,,FF::FF:F             |
| **[XB](glossary.html#xb_auxiliary_tag)**  | The probability that sample barcode decoding is incorrect.     | 2.27479e-06                   |
| **[CB](glossary.html#cb_auxiliary_tag)**  | Cellular identifier.                                           | ACTGCATA                      |
| **[CR](glossary.html#cr_auxiliary_tag)**  | Raw uncorrected cellular barcode sequence.                     | ACTGCATT                      |
| **[CY](glossary.html#cr_auxiliary_tag)**  | Phred quality of the cellular barcode sequence in the CR tag.  | ,,FF::FF                      |
| **[XC](glossary.html#xc_auxiliary_tag)**  | The probability that Cellular barcode decoding is incorrect.   | 2.27479e-06                   |
| **[MI](glossary.html#mi_auxiliary_tag)**  | Molecular Identifier.                                          |                               |
| **[RX](glossary.html#rx_auxiliary_tag)**  | Molecular barcode sequence, either corrected or uncorrected.   |                               |
| **[QX](glossary.html#qx_auxiliary_tag)**  | Phred quality of the molecular barcode sequence in the RX tag. |                               |
| **[OX](glossary.html#ox_auxiliary_tag)**  | Raw uncorrected molecular barcode sequence.                    |                               |
| **[BZ](glossary.html#bz_auxiliary_tag)**  | Phred quality of the molecular barcode sequence in the OX tag. |                               |
| **[XM](glossary.html#xm_auxiliary_tag)**  | The probability that molecular barcode decoding is incorrect.  |                               |

---

## Key Concepts
