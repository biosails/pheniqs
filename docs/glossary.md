---
layout: default
title: "Glossary"
permalink: /glossary
id: glossary
---

<a name="zero_based_coordinate" />
1 based coordinate system
: A coordinate system where the first base of a sequence is one. In this coordinate system both start and end coordinates are inclusive. The total length of a token is given by `end - start + 1`.

<a name="one_based_coordinate" />
0 based coordinate system
: A coordinate system where the first base of a sequence is zero. In this coordinate system, the start coordinate in inclusive and the end coordinate is exclusive. The total length of a token is given by `end - start`.

<a name="template_sequence" />Template sequence
: A DNA/RNA sequence part of which is sequenced on a sequencing machine or assembled from raw sequences.

<a name="technical_sequence" />Technical sequence
: A synthetic nucleotide sequence often ligated to a [template sequence](#template_sequence) for the purpose of classification, quantification or some other quality control need. Some common examples of a technical sequence are a sequencing adapter, [multiplex barcode](#multiplex_barcode) or a [molecular barcode](#molecular_barcode).

<a name="read" />Sequence Read
: A raw sequence that comes off a sequencing machine. A read may consist of multiple ordered [segments](#segment).

<a name="segment" />Sequence Read Segment
: A contiguous [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) encoded sequence optionally accompanied by the corresponding [Phred](https://en.wikipedia.org/wiki/Phred_quality_score) encoded quality scores and [auxiliary meta data tags](http://samtools.github.io/hts-specs/SAMtags.pdf). A segment can contain a [template sequence](#template_sequence), [technical sequence](#technical_sequence) or a combination of the two.

<a name="input_segment" />Input segment
: An enumerated [segment](#segment) of the input [read](#read). Input segment indexes are [zero based](#zero_based_coordinate).

<a name="output_segment" />Output segment
: An enumerated [segment](#segment) of the output [read](#read). Output segment indexes are [zero based](#zero_based_coordinate).

<a name="multiplex_barcode" />Multiplex barcode
: A [technical sequence](#technical_sequence) used to classify a [read](#read) into read groups. A multiplex barcode can have multiple [segments](#segment) and a match to all segments is required by the decoder to declare a successful classification. The sequence and quality of the concatenated multiplex barcode are written in SAM records to the standardized [BC](#bc_auxiliary_tag) and [QT](#qt_auxiliary_tag) auxiliary tags. When decoding with [Phred-adjusted maximum likelihood](#phred_adjusted_maximum_likelihood_decoding) Pheniqs also writes the decoding error probability to the [XB](#dq_auxiliary_tag) tag. Multiplex barcode segment indexes are [zero based](#zero_based_coordinate).

<a name="molecular_barcode" />Molecular barcode
: Also known as a UMI (Unique Molecular Identifier) is a [technical sequence](#technical_sequence) used to identify unique DNA fragments. Duplicate read elimination is often used to assist in quantification and to improve read confidence. The sequence and quality of the molecular barcode are written in SAM records to the community proposed [RX](#rx_auxiliary_tag) and [QX](#qx_auxiliary_tag) auxiliary tags.

<a name="split_file_layout" />Split file layout
: Read [segments](#segment) are split into separate files by [read group](#read_group) and segment index. All segments in a file have the same read group. A file contains only a single segment from each read. This is a common layout for [FASTQ](#fastq) files.

<a name="interleaved_file_layout" />Interleaved file layout
: Read [segments](#segment) are split into separate files by [read group](#read_group). All segments in a file have the same read group. A file contains multiple segments from each read. Read segments are consecutive and in order which implies the GO property of the [@HD header tag](#hd_header_tag) is set to `query`. This is a common layout for SAM based formats but is also sometimes used with [FASTQ](#fastq).

<a name="combined_file_layout" />Combined file layout
: Read [segments](#segment) from multiple [read groups](#read_group) are written to the same file. A file contains multiple segments from each read. Read segments are consecutive and in order which implies the GO property of the [@HD header tag](#hd_header_tag) is set to `query`. This layout can be more efficient for SAM based formats. Although writing combined output to [FASTQ](#fastq) is supported, the lack of standardized meta data encoding in FASTQ makes this impractical in real world scenarios.

<a name="read_group" />Read Group
: For detailed semantics of the various read group fields you should consult the [sequence alignment map format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="fi_auxiliary_tag" />FI auxiliary tag
: The [one based](#one_based_coordinate) index of a [segment](#segment) in the [read](#read).

<a name="tc_auxiliary_tag" />TC auxiliary tag
: The number of [segments](#segment) in the [read](#read).

<a name="rg_auxiliary_tag" />RG auxiliary tag
: Read group identifier. The value must be consistent with the ID property of an [@RG](#rg_header_tag) header tag.

<a name="bc_auxiliary_tag" />BC auxiliary tag
: Raw uncorrected multiplex barcode sequence with any quality scores stored in the [QT](#qt_auxiliary_tag) tag. The BC tag should match the QT tag in length. In the case of multiple multiplex barcode segments, all segments are concatenated with a hyphen (‘-’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="qt_auxiliary_tag" />QT auxiliary tag
: [Phred](#sanger_format) quality of the sample barcode sequence in the [BC](#bc_auxiliary_tag) tag. Phred score is + 33 encoded. In the case of multiple multiplex barcode segments, all segments are concatenated with a space (‘ ’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="xb_auxiliary_tag" />XB auxiliary tag
: The probability that the decoding of the multiplexing barcode stored in [BC](#bc_auxiliary_tag) tag was incorrect.

<a name="rx_auxiliary_tag" />RX auxiliary tag
: Sequence bases from the [unique molecular identifier](#molecular_barcode) with any quality scores stored in the [QX](#qx_auxiliary_tag) tag. The RX tag should match the QX tag in length. These could be either corrected or uncorrected. Unlike [MI](#mi_auxiliary_tag) tag, the value may be non-unique in the file. In the case of multiple molecular barcode segments, all segments are concatenated with a hyphen (‘-’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf). If the bases represent corrected bases, the original sequence can be stored in [OX](#ox_auxiliary_tag) tag.

<a name="qx_auxiliary_tag" />QX auxiliary tag
: [Phred](#sanger_format) quality of the [unique molecular identifier](#molecular_barcode) sequence in the [RX](#rx_auxiliary_tag) tag. Phred score + 33 encoded. The qualities here may have been corrected with raw bases and qualities stored in [OX](#ox_auxiliary_tag) tag and [BZ](#bz_auxiliary_tag) tag respectively. In the case of multiple molecular barcode segments, all segments are concatenated with a space (‘ ’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf). If the qualities represent corrected values, the original values can be stored in [BZ](#bz_auxiliary_tag) tag.

<a name="ox_auxiliary_tag" />OX auxiliary tag
: Raw uncorrected [unique molecular identifier](#molecular_barcode) bases, with any quality scores stored in the [BZ](#bz_auxiliary_tag) tag. In the case of multiple molecular barcode segments, all  segments are concatenated with a space (‘ ’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="bz_auxiliary_tag" />BZ auxiliary tag
: [Phred](#sanger_format) quality of the uncorrected unique molecular identifier sequence in the [OX](#ox_auxiliary_tag) tag. Phred score + 33 encoded. In the case of multiple molecular barcode segments, all segments are concatenated with a space (‘ ’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="mi_auxiliary_tag" />MI auxiliary tag
: Molecular Identifier. A unique ID within the SAM file for the source molecule from which this read is derived. All reads with the same MI tag represent the group of reads derived from the same source molecule.

<a name="xm_auxiliary_tag" />XM auxiliary tag
: The probability that the decoding of the molecular barcode stored in [RX](#rx_auxiliary_tag) tag was incorrect.

<a name="cb_auxiliary_tag" />CB auxiliary tag
: Unique cell identifier. Pheniqs populates this tag with the corrected sequence bases of the cellular barcode.

<a name="cr_auxiliary_tag" />CR auxiliary tag
: Raw uncorrected cellular identifier bases, with any quality scores stored in the [CY](#cy_auxiliary_tag) tag. In the case of multiple cellular barcode segments, all segments are concatenated with a hyphen (‘-’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="cy_auxiliary_tag" />CY auxiliary tag
: [Phred](#sanger_format) quality of the uncorrected cellular identifier sequence in the [CY](#cy_auxiliary_tag) tag. Phred score + 33 encoded. In the case of multiple cellular barcode segments, all segments are concatenated with a space (‘ ’) separator, following the sequence alignment map format specification [recommendation](https://samtools.github.io/hts-specs/SAMv1.pdf).

<a name="xc_auxiliary_tag" />XC auxiliary tag
: The probability that the decoding of the cellular barcode stored in [CB](#cb_auxiliary_tag) tag was incorrect.

<a name="xe_auxiliary_tag" />XE auxiliary tag
: The expected number of base call errors in the [segment](#segment) as computed from the quality scores.

<a name="lb_auxiliary_tag" />LB auxiliary tag
: Library name.

<a name="sm_auxiliary_tag" />SM auxiliary tag
: Sample name.

<a name="pg_auxiliary_tag" />PG auxiliary tag
: Program record identifier. The value must be consistent with the ID property of an [@PG](#pg_header_tag).

<a name="pu_auxiliary_tag" />PU auxiliary tag
: Platform unit. The value must be consistent with the PU property of an [@RG](#rg_header_tag) header tag.

<a name="co_auxiliary_tag" />CO auxiliary tag
: Free-text comment.

<a name="hd_header_tag" />@HD header tag
: Defines the format version, record sort order and grouping in a SAM header.

<a name="rg_header_tag" />@RG header tag
:   Defines a read group in a SAM header.

<a name="rg_id_header_tag" />@RG:ID header tag
:   Read group identifier. Each @RG line must have a unique ID. The value of ID is used in the RG tags of alignment records. Must be unique among all read groups in header section.

<a name="rg_bc_header_tag" />@RG:BC header tag
:   Barcode sequence identifying the sample or library. This value is the expected barcode bases as read by the sequencing machine in the absence of errors.

<a name="rg_cn_header_tag" />@RG:CN header tag
:   Name of sequencing center producing the read.

<a name="rg_ds_header_tag" />@RG:DS header tag
:   Description. UTF-8 encoding may be used.

<a name="rg_dt_header_tag" />@RG:DT header tag
:   Date the run was produced (ISO8601 date or date/time).

<a name="rg_lb_header_tag" />@RG:LB header tag
:   Library associated with the read group.

<a name="rg_pg_header_tag" />@RG:PG header tag
:   Programs used for processing the read group.

<a name="rg_pi_header_tag" />@RG:PI header tag
:   Predicted median insert size.

<a name="rg_pl_header_tag" />@RG:PL header tag
:   Platform/technology used to produce the reads. Valid values: *CAPILLARY*, *DNBSEQ* (MGI/BGI), *HELICOS*, *ILLUMINA*, *IONTORRENT*, *LS454*, *ONT* (Oxford Nanopore), *PACBIO* (Pacific Biotechnology), and *SOLID*. Pheniqs uses this tag to interpret metadata from the read ID in FASTQ formatted reads.

<a name="rg_pm_header_tag" />@RG:PM header tag
:   Platform model. Free-form text providing further details of the platform/technology used.

<a name="rg_pu_header_tag" />@RG:PU header tag
:   Platform unit unique identifier.

<a name="rg_sm_header_tag" />@RG:SM header tag
:   Platform unit unique identifier.

<a name="pg_header_tag" />@PG header tag
: Defines a program in a SAM header.

<a name="phred_adjusted_maximum_likelihood_decoding" />Phred adjusted maximum likelihood decoding
: A maximum likelihood decoder that directly estimates the decoding error probability from the base calling error probabilities provided by the sequencing platform. Abbreviated **PAMLD**.

<a name="minimum_distance_decoding" />Minimum distance decoding
: A minimum distance decoder, abbreviated **MDD**, consults the edit distance between the expected and observed sequence. It does not take base calling error probabilities into consideration. For more information see [minimum distance decoding on wikipedia](https://en.wikipedia.org/wiki/Decoding_methods#Minimum_distance_decoding). For a good review of the application of minimum distance decoding to sequencing see
*[Short Barcodes for Next Generation Sequencing](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082933) by Mir, K. et al*.

<a name="htslib" />HTSlib
: [HTSlib](http://www.htslib.org) is an Open Source implementation of a unified multi-threaded C library for accessing high-throughput sequencing data encoded according to the [sequence alignment map format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) in either the SAM, BAM or CRAM format. It also provides efficient and standardized means of manipulating [FASTQ](#fastq) files. SAM is a human readable TAB-delimited text format superior to FASTQ in the sense that it supports standardized meta data annotations and is easier to parse. BAM is a BGZF compressed binary encoding of SAM data that provides efficient random access for indexed queries. [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf) is a newer, more advanced, binary encoding that is fully compatible with BAM and offers significantly better lossless compression, a lossy compression scheme and improved IO but does not support encoding degenerate [IUPAC](https://en.wikipedia.org/wiki/Nucleic_acid_notation) nucleotides. HTSlib is licensed under the [MIT license](https://opensource.org/licenses/MIT).

<a name="fastq" />FASTQ
: a text-based file format for storing both a biological sequence and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity. [SAM](#htslib) encoded files have several major advantages over FASTQ files like standardized metadata tags, better compression, efficient random access and easier parsing and are slowly replacing FASTQ files. See the [FASTQ wikipedia page](https://en.wikipedia.org/wiki/FASTQ_format) for more details about the FASTQ file format.

<a name="sanger_format" />Sanger Format
: A Phred quality encoding that encodes quality scores from 0 to 93 using ASCII 33 to 126. Since encoding to ASCII involves adding 33 to the Phred value the encoding offset is said to be **33**. Current Illumina platforms directly produce FASTQ in Sanger format and it is the only quality encoding allowed in [SAM](#htslib) records. See the [FASTQ wikipedia page for more details](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).

<a name="qc_fail" />QC Fail Flag
: Some reads are marked as failing quality control. This is [signaled](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) on the comment portion of the read identifier in [FASTQ](glossary#fastq) files or the **512** flag on a SAM record flag. Reads that pass quality control are called **PF** for **pass-filter** (SAM flag **ON**), while those that fail are often referred to as **QC fail** reads (SAM flag **OFF**).

<a name="relative_path" />Relative path
: A relative path is one that **does not** begin with `/`. Relative paths must be resolved against a base directory path. By default Pheniqs resolves those against work directory where Pheniqs was executed unless `base input url` or `base output url` was specified in which case input and output paths are resolved against those base directories.

<a name="absolute_path" />Absolute path
: An absolute path is one that begin with `/`. Absolute paths ignore `base input url` or `base output url`.

<a name="standard_stream" />Standard stream
: POSIX systems use standard streams, stdin, stdout and stderr to transfer data between applications without writing it to a file on a hard drive. This allows to pipe the output from one utility into the input of the next in the pipeline without the need for additional storage. Avoiding writing to the disk can also significantly accelerate some pipelines since drives are often orders of magnitude slower than computer memory.

<a name="static_linking" />Static linking
: A [statically linked](https://en.wikipedia.org/wiki/Static_library) build copies all or most of the library dependencies into the binary executable, producing a standalone executable that can be moved around the system and to compatible systems. This method, although yielding a larger executable, resolves conflicts with different version of the dependencies that might already be present on the system. Building a portable binary from scratch using the built-in package manager makes any Pheniqs revision easily available on cluster environments where the user rarely has elevated permissions required to perform a standard build.

<a name="feed_resolution" />Feed resolution
: The number of consecutive read segments in the feed that have the same read identifier. In most interleaved files all segments of a read are adjacent and in the same segment order.

<!--
<a name="" />
TERM
: DEFINITION
 -->
