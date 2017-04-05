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
        <li><a                  href="/pheniqs/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/manual.html">Documentation</a></li>
        <li><a class="active"   href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="https://github.com/GunsalusPiano/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Glossary
{:.page-title}

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
: A [technical sequence](#technical_sequence) used to classify a [read](#read) into read groups. A multiplex barcode can have multiple [segments](#segment) and a match to all segments is required by the decoder to declare a successful classification. The sequence and quality of the concatenated multiplex barcode are written in SAM records to the standardized [BC](#bc_auxiliary_tag) and [QT](#qt_auxiliary_tag) auxiliary tags. When decoding with [Phred-adjusted maximum likelihood](#phred_adjusted_maximum_likelihood_decoding) Pheniqs also writes the decoding error probability to the [DQ](#dq_auxiliary_tag) tag. Multiplex barcode segment indexes are [zero based](#zero_based_coordinate).

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
: Multiplex barcode sequence, with quality scores stored in the [QT](#qt_auxiliary_tag) tag.

<a name="qt_auxiliary_tag" />QT auxiliary tag
: [Phred](#sanger_format) encoded quality of the multiplex barcode sequence in the [BC](#bc_auxiliary_tag) tag.

<a name="rx_auxiliary_tag" />RX auxiliary tag
: [Molecular barcode](#molecular_barcode) sequence, with quality scores stored in the [QX](#qx_auxiliary_tag) tag.

<a name="qx_auxiliary_tag" />QX auxiliary tag
: [Phred](#sanger_format) encoded quality of the [molecular barcode](#molecular_barcode) sequence in the [RX](#rx_auxiliary_tag) tag.

<a name="dq_auxiliary_tag" />DQ auxiliary tag
: A proposed auxiliary tag of type float to store the [multiplex barcode](#multiplex_barcode) decoding error probability.

<a name="ee_auxiliary_tag" />EE auxiliary tag
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

<!-- 
<a name="" />
TERM
: DEFINITION
 -->
