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
        <li><a class="active"   href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Build</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Vignettes
{:.page-title}

* placeholder
{:toc}

Pheniqs can be configured to handle any arbitrary configuration of biological and technical sequences such as barcoded libraries, cellular indexes, and UMIs, for both bulk and single-cell experiments. \textit{Transform patterns} define the number, location, and type of barcodes to be extracted, which will vary depending on the specific experimental design. Different types of sequences are extracted using \emph{tokens}:

![transform patterns](transform_patterns.png)

See the [**Overview**](overview.html) page for more information on common experimental designs, read anatomy, tokenization, and transform patterns. The [**Tokenization**](transform.html) page provides a high-level view of how transform patterns are constructed and the output they provide.

Examples of how to configure Pheniqs for a handful of published experimental designs are provided here to help users get started with their own applications.

# Background

To configure Pheniqs for any particular workflow, it helps to have some idea about next-gen sequencing technology, Illumina's bcl2fastq softare, common file formats, and some general terminology.

<!-- Right now this is a laundry list that needs to be organized. -->

## Illumina Sequencing

Illumina offers several sequencing platforms for short-read sequencing that have different throughput capacity and somewhat different chemistry.

+ An Introduction to Next-Generation Sequencing Technology [PDF](https://www.illumina.com/documents/products/illumina_sequencing_introduction.pdf)


### Sequence file formats

Pheniqs can accept as input and write as output both FASTQ and SAM/BAM/CRAM sequence formats. In addition, Pheniqs can perform on-the-fly conversion between these sequence formats.

#### FASTQ

+ [FASTQ  format](https://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm)
+ [FASTQ files explained](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)

#### Sequence Alignment/Map (SAM)

+ SAM specification [PDF](https://samtools.github.io/hts-specs/SAMv1.pdf)
+ Auxilliary tags [PDF](https://samtools.github.io/hts-specs/SAMtags.pdf)


### bcl2fastq

bcl2fastq is used to convert raw BCL files produced by Illumina sequencers into FASTQ files for input to Pheniqs.

+ bcl2fastq2 Conversion Software v2.20 [PDF](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf)

### Sample Sheets

Metadata about any run is contained in Illumina **sample sheets**. Pheniqs is packaged with [Python helper tools](pyapi.html) to assist users with I/O management. The [IlluminaAPI](illumina_python_api.html) can automatically generate an initial configuration file from any Illumina sample sheet.

+ [More information about Illumina sample sheets]()

Barcode design: To optimize demultiplexing results, choose index adapters that optimize color balance when performing library prep.




### Adapter Sequences

Sequence libraries are prepared with flanking **adapter sequences** that are used as anchors by the sequencing machine during sequencing.

+ [Indexed Sequencing Overview Guide](https://support.illumina.com/downloads/indexed-sequencing-overview-15057455.html) - [PDF](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-05.pdf) (updated April 2, 2019)
+ [Illumina Adapter Sequences]](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html) - [PDF](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-11.pdf) (updated June 13, 2019)
+ [IDT adapters](https://www.idtdna.com/pages/products/next-generation-sequencing/adapters)
  - Article: Unique, dual-indexed sequencing adapters with UMIs effectively eliminate index cross-talk and significantly improve sensitivity of massively parallel sequencing. MacConaill LE, Burns RT, et al. [ _BMC Genomics_, 19 : 30. (2018) ](https://www.idtdna.com/pages/products/next-generation-sequencing/adapters)

#### Adapter Trimming and UMI Removal

Depending on settings, the bcl2fastq2 Conversion Software can trim adapter sequences and remove Unique Molecular Identifier (UMI) bases from reads. _**UMI removal by bcl2fastq is not recommended since Pheniqs can extract and preserve all barcodes in SAM format.**_

+ Adapter trimming — The software determines whether a read extends past the DNA insert and into the sequencing adapter. An approximate string matching algorithm identifies all or part of the adapter sequence and treats inserts and deletions (indels) as one mismatch. Base calls matching the adapter sequence and beyond are masked or removed from the FASTQ file.
+ UMI removal — UMIs are random k-mers attached to the genomic DNA (gDNA) before polymerase chain reaction (PCR) amplification. After the UMI is amplified with amplicons, the software can retrieve the bases and include them in the read name in the FASTQ files. When the TrimUMI sample sheet setting is active, the software can also remove the bases from the reads.


## Some Terminology

Based on: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups

In practice, a **read group** refers to a set of reads that were generated from a single run of a sequencing instrument. Read groups are identified in the SAM/BAM/CRAM file by a number of tags that are defined in the official SAM specification.

In the simple case where a single library preparation derived from a single biological sample was run on a single lane of a flowcell, all the reads from that lane run belong to the same read group.

When multiplexing is involved, then each subset of reads originating from a separate library run on the same lane will constitute a separate read group.

Some analysis tools will fail if the read groups are not specified in the BAM file's header section of the BAM file. The following command shows the read group information: `samtools view -H sample.bam | grep @RG`

Additional information about libraries and samples are also specified by specific SAM tags. The basic set of tags associate with read groups are as follows:

#### ID = Read group identifier
This tag identifies which read group each read belongs to, so each read group's ID must be unique. It is referenced both in the read group definition line in the file header (starting with @RG) and in the RG:Z tag for each read record. In Illumina data, read group IDs are composed using the flowcell + lane name and number, making them a globally unique identifier across all sequencing data in the world.

Note that some tools have the ability to modify IDs when merging SAM files in order to avoid collisions. ID is the lowest denominator that differentiates factors contributing to technical batch effects: therefore, a read group is effectively treated as a separate run of the instrument in data processing steps such as base quality score recalibration, since they are assumed to share the same error model.

#### PU = Platform Unit
The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. The {FLOWCELL_BARCODE} refers to the unique identifier for a particular flow cell. The {LANE} indicates the lane of the flow cell and the {SAMPLE_BARCODE} is a sample/library-specific identifier.

Note that this field can only be added after demultiplexing.

SM = Sample
The name of the sample sequenced in this read group.
In the case of multiplexing, in which multiple samples are run together in one lane, the SM tag defines the individual sample names.
Otherwise there will be one SM tag per lane.
PL = Platform/technology used to produce the read
This constitutes the only way to know what sequencing technology was used to generate the sequencing data. Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.
LB = DNA preparation library identifier
All samples derived from the same library prep should share the same LB tag. The LB field helps determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes.
The above document illustrates a vignette for a case where samples from three individuals were used to prepare two different libraries that were each run across two lanes of an Illumina HiSeq, generating 3 x 2 x 2 = 12 read groups that are written to 12 BAM files (one read group / lane / BAM file). There is a hierarchical relationship between read groups (unique for each lane) to libraries (sequenced on two lanes) and samples (across four lanes, two lanes for each library). Each read group is assigned a unique ID and carries additional metadata:
Example: @RG ID:FLOWCELL1.LANE1 PL:Illumina LB:LIB1 SM:SAMPLE1 PI:200
ID: Flowcell#.Lane# (unique per lane)
PL: Illumina (common for all files)
LB: Unique per library (here, shared between 2 lanes each)
SM: Unique per sample; shared between all libraries for each sample (here, 2 libraries x 2 lanes per sample = 4 read groups)
PI: 200 or 400 (fragment lengths)

Note from GATK documentation: GATK expects all read groups appearing in the read data to be specified in the file header, and will fail with an error if it does not find that information (whether there is no read group information in the file, or a subset of reads do not have read groups).

For BAM files that lack required fields or do not differentiate pertinent factors within the fields, Picard's AddOrReplaceReadGroups can be used to add or appropriately rename the read group fields as outlined here.

---

# Indexing Strategies

**Barcodes** allow the disambiguation of distinct experimental datasets such as samples, libraries, and individual cells. Barcodes are composed of one or more index sequences that can be added in different locations in a PCR product and can be joined in a variety of ways to form combinatorial barcodes. Barcoding may be performed using multiplex or inline barcodes, or a mixture of these.

**Multiplex index sets** are provided by several vendors to enable single- or dual-indexing of samples. On the Illumina platform, these are sequenced separately during sequencing and are output in the Index 1 and Index 2 reads. These optional tags require additional rounds of priming, are emitted as separate read segments during sequencing, and are recognized and treated as barcode indexes by standard analysis pipelines.

An advantage to using multiplex barcodes is that they don't take up space within the sample reads, allowing more biological sequence to be recorded. This is especially helpful for short read runs and for applications where every base is precious, such as de novo assembly of novel genomes.

**Inline barcodes**, instead, are either ligated to the sample DNA or are encoded in primers used during the RT (reverse transcription) or amplification steps. They take up some portion of the full-length read and need to be extracted during a post-processing step.

An advantage of using inline barcodes is that they accommodate a wide variety of experimental designs and any number of indexes can be combined to form compound barcodes comprising two or more individual index tags. Inline barcodes are convenient for applications that work with well annotated genomes, such as RNA-seq, ChIP-seq, ATAC-seq, SNP calling, etc.

Due to practical considerations of library preparation, inline barcodes are usually located at the beginning of Read 1 or Read 2. Because the sequencer performs color calibration at the beginning of the reads, inline barcodes should be as base diverse as possible.


# Vignettes




Pheniqs can accommodate virtually any experimental design due to its flexible syntax for parsing read segments. Some common designs for the Illumina platform are illustrated here:

![read anatomy](/pheniqs/assets/img/diagram8.png)
<a name="illumina_python_api" />
[Prior estimated Illumina with the python API](illumina_python_api.html)
: This vignette will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the Pheniqs python API. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors and demultiplex the run. It loosely applies to almost every standard sample multiplex Illumina run.

<a name="standard_illumina" />
[Standard Illumina sample demultiplexing](illumina.html)
: This vignette will walk you through writing configuration files for demultiplexing a standard Illumina high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
[Fluidigm with a sample and a cellular tag](fluidigm.html)
: This vignette will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

# Read anatomy

Illumina sequencing platforms typically produce four different sequence elements: two Index sequences, referred by Illumina as the **i5** and **i7** barcodes, and two Insert sequences, referred by Illumina as **read 1** and **read 2**. Collectively, these are referred to as read segments. For example, consider a [standard paired-end, dual index library design](illumina.html):

![read anatomy](/pheniqs/assets/img/diagram1.png)
<!-- <img src="/pheniqs/assets/img/diagram1.png" width="75%" class="figure_medium" /> -->

The read segments for this standard design thus comprise two technical sequences (referred by Illumina as I1, I2) and two biological sequences (referred by Illumina as R1, R2):

<img src="/pheniqs/assets/img/diagram2.png" width="75%" class="figure_medium" />

# Sequence Classification

The combination of the barcodes contained in the **I1** and **I2** index positions specifies the sample library. With standard dual indexing, up to 96 distinct sample libraries can be pooled and run together in a single sequencing lane.

To identify which biological sequences belong to which library, the sequences belonging each one need to be separated from each other. This process of deconvolving libraries is called demultiplexing and is done by classifying each of the sequences using the barcode indexes:

![read anatomy](/pheniqs/assets/img/diagram5.png)

>**Note** The i5 adaptor sequences specified in sample sheets will be reverse complemented for platforms that read the I2 index on the bottom strand. The i7 sequences in sample sheets are always reverse complemented relative to the original adaptor sequences since they are read from the top strand
{: .example}

# Tokenization

Pheniqs uses [tokens](manual.html#tokenization) to reference and extract information from different read segments by specifying where to look for different classes of sequence elements (i.e. barcodes, biological sequences). Each element of interest is defined by an offset relative to the beginning of a given read segment (in this example I1, I2, R1, R2) and a length. It is important to note that Pheniqs uses [zero based](glossary.html#zero_based_coordinate) indexing, so the first read to come off the machine will be referred to as Segment 0, and so on:

![read anatomy](/pheniqs/assets/img/diagram7.png)

For a standard paired-end, dual indexed Illumina run, the sample barcodes usually comprise the full I1 and I2 read segments. Because Illumina sequencing is calibrated in relation to the previously sequenced base, those segments are sometimes sequenced one nucleotide longer than necessary to ensure good quality on the last nucleotide. The biological sequences start at the first position of R1 and R2 and extend for the full number of cycles run (typically 75, 100, or 150 nucleotides).

For this design, the barcode tokens begin at position 0 in I1 and I2 and extend for 10 bases. The tokens for biological sequences begin at position 0 of Read 1 and Read 2 and extend for the full span of those read segments.

# Input / Output

Pheniqs can manipulate [SAM, BAM and CRAM](glossary.html#htslib) files as well as uncompressed and gzip compressed [FASTQ](glossary.html#fastq). Configuration and reports are [JSON](https://en.wikipedia.org/wiki/JSON) encoded for easy integration.

All of these can be specified with standard templates and may be overridden by additional directives within the configuration file.

After sequence classification, the various tokenized sequences extracted are written along with their respective confidence scores to specific [auxiliary tags](https://samtools.github.io/hts-specs/SAMtags.pdf) within the [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format, as outlined below.

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
