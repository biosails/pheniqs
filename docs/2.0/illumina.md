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
        <li><a                  href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Build</a></li>
        <li><a class="active"   href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/best_practices.html">Best Practice</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Standard Illumina sample decoding
{:.page-title}

This tutorial demonstrates how to manually prepare configuration files for decoding sample barcodes for a standard Illumina sequencing run. Pheniqs includes a Python API that helps users create configuration files automatically. For a tutorial that uses the Python API to generate the configuration files for this example, see the [Standard Illumina sample decoding with the python API](illumina_python_api.html).

In this example the run has paired-end dual-index samples multiplexed using the standard Illumina i5 and i7 index protocol. The read is made of 4 segments: 2 biological sequences (cDNA, genomic DNA, etc.) read from both ends of the insert fragment, and 2 technical sequences containing the i5 and i7 indices. If the results are written to SAM, BAM or CRAM the multiplex barcode and its quality scores are written to the [BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag) tags, and the decoding error probabilities are written to the [XB](glossary.html#xb_auxiliary_tag) tag.

## Input Read Layout

![Illumina paired-end dual-index sequencing](/pheniqs/assets/img/Illumina_paired-end_dual-index.png)

Base calling with bc2fastq will produce 4 files per lane:
- `L001_R1_001.fastq.gz`: Read 1, starting from the beginning of the insert fragment ("top" strand).
- `L001_I1_001.fastq.gz`: Index 1, the i7 index.
- `L001_I2_001.fastq.gz`: Index 2, the i5 index.
- `L001_R2_001.fastq.gz`: Read 2, starting from the other end of the insert fragment ("bottom" strand); note that this sequence is reverse complemented relative to the Read 1 sequence.

>```json
"input": [
    "Lane1_S1_L001_R1_001.fastq.gz",
    "Lane1_S1_L001_I1_001.fastq.gz",
    "Lane1_S1_L001_I2_001.fastq.gz",
    "Lane1_S1_L001_R2_001.fastq.gz"
],
```
>**declaring input read segments** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling.
{: .example}

To emit the two ends of the insert region as two segments of the output read, we declare the global transform directives
>```json
"transform": { "token": [ "0::", "3::" ] }
```
>**declaring output read segments** Only the segments coming from the first and the fourth file are biological sequences and should be included in the output.
{: .example}

To classify the reads by the i5 and i7 indices, we declare a `codec`, which is the list of possible barcode sequences, and a `transform` that tells Pheniqs which read segment(s) and coordinates correspond to the barcode sequence(s).
>```json
"multiplex": {
  "comment": "Sample barcodes are 96 unique combinations of i5 (10nt) + i7 (10nt) sequences.",
  "transform": { "token": [ "1::10", "2::10" ] },
  "codec": {
    "@GAACTGAGCGCGCTCCACGA": { "barcode": [ "GAACTGAGCG", "CGCTCCACGA" ] },
    "@GACGAGATTAAGGATAATGT": { "barcode": [ "GACGAGATTA", "AGGATAATGT" ] }
  }
}
```
>**declaring sample demultiplexing** The standard Illumina dual-index protocol allows up to 96 unique dual 10 base barcodes in the i5 and i7 region, but for the sake of brevity we only show 2 here. Comments are allowed within any dictionary element in the configuration file and are ignored by Pheniqs. Alternatively, a codec may also be inherited from a base decoder rather than being enumerated explicitly here (see **inheritance**).
{: .example}

To discard reads that failed the internal Illumina sequencer noise filter, we instruct Pheniqs to filter incoming *qc fail* reads
>```json
"filter incoming qc fail": true
```

Putting things together, we can generate a configuration for demultiplexing the reads into three output files using the *PAMLD* decoder. Here, we declare that we expect **5%** of reads to be foreign DNA, which should not classify to any of the barcodes (for instance spiked-in PhiX sequences). We specify the noise prior with the `"noise": 0.05` directive. To mark as *qc fail* those reads with a posterior probability of correct barcode decoding that is below **0.95**, we specify `"confidence threshold": 0.95`.

>```json
{
    "CN": "CGSB",
    "PL": "ILLUMINA",
    "PM": "HiSeq",
    "base input url": "~/CBJLFACXX/raw",
    "base output url": "~/CBJLFACXX/sample",
    "filter incoming qc fail": true,
    "flowcell id": "CBJLFACXX",
    "input": [
        "Lane1_S1_L001_R1_001.fastq.gz",
        "Lane1_S1_L001_I1_001.fastq.gz",
        "Lane1_S1_L001_I2_001.fastq.gz",
        "Lane1_S1_L001_R2_001.fastq.gz"
    ],
    "transform": { "token": [ "0::", "3::" ] },
    "multiplex": {
      "transform": { "token": [ "1::10", "2::10" ] },
      "algorithm": "pamld",
      "noise": 0.05,
      "confidence threshold": 0.95,
      "codec": {
          "@GAACTGAGCGCGCTCCACGA": {
              "SM": "c57bl6 mouse",
              "LB": "c57bl6 fetal liver",
              "barcode": [ "GAACTGAGCG", "CGCTCCACGA" ]
          },
          "@GACGAGATTAAGGATAATGT": {
              "SM": "c57bl6 mouse",
              "LB": "c57bl6 adult bone marrow",
              "barcode": [ "GACGAGATTA", "AGGATAATGT" ]
          }
      }
    }
}
```
>**Single index Paired end Illumina protocol** Classifying the 2 barcodes using the 4 fastq files produced by bcl2fastq.
{: .example}

Before we proceed we validate the configuration with Pheniqs:

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json --validate
```

The output is a readable description of all the explicit and implicit parameters after applying defaults. You can check how Pheniqs detects the input format, compression and layout as well as the output you can expect. In The *Output transform* section is a verbal description of how the output read segments will be assembled from the input. similarly *Transform* in the *Mutliplex decoding* section desribes how the segment that will be matched against the barcodes is assembled. the You can also see how each of the read groups will be tagged and the prior probability PAMLD will assume for each barcode.

    Environment

        Base input URL                              /home/lg/CBJLFACXX/raw
        Base output URL                             /home/lg/CBJLFACXX/raw
        Platform                                    ILLUMINA
        Quality tracking                            disabled
        Filter incoming QC failed reads             enabled
        Filter outgoing QC failed reads             disabled
        Input Phred offset                          33
        Output Phred offset                         33
        Leading segment index                       0
        Default output format                       sam
        Default output compression                  none
        Default output compression level            5
        Feed buffer capacity                        2048
        Threads                                     8
        Decoding threads                            1
        HTSLib threads                              8

    Input

        Input segment cardinality                   4

        Input segment No.0 : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz
        Input segment No.1 : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz
        Input segment No.2 : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz
        Input segment No.3 : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz

        Input feed No.0
            Type : fastq
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz

        Input feed No.1
            Type : fastq
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz

        Input feed No.2
            Type : fastq
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz

        Input feed No.3
            Type : fastq
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/CBJLFACXX/raw/Lane1_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz

    Output transform

        Output segment cardinality                  2

        Token No.0
            Length        variable
            Pattern       0::
            Description   cycles 0 to end of input segment 0

        Token No.1
            Length        variable
            Pattern       3::
            Description   cycles 0 to end of input segment 3

        Assembly instruction
            Append token 0 of input segment 0 to output segment 0
            Append token 1 of input segment 3 to output segment 1

    Mutliplex decoding

        Decoding algorithm                          pamld
        Shannon bound                               2 3
        Noise                                       0.05
        Confidence threshold                        0.95
        Segment cardinality                         2
        Nucleotide cardinality                      20
        Barcode segment length                      10 10

        Transform

            Token No.0
                Length        10
                Pattern       1::10
                Description   cycles 0 to 10 of input segment 1

            Token No.1
                Length        10
                Pattern       2::10
                Description   cycles 0 to 10 of input segment 2

            Assembly instruction
                Append token 0 of input segment 1 to output segment 0
                Append token 1 of input segment 2 to output segment 1


        Barcode undetermined
            ID : CBJLFACXX:undetermined
            PU : CBJLFACXX:undetermined
            PL : ILLUMINA
            PM : HiSeq
            CN : CGSB
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Barcode @GAACTGAGCGCGCTCCACGA
            ID : CBJLFACXX:GAACTGAGCGCGCTCCACGA
            PU : CBJLFACXX:GAACTGAGCGCGCTCCACGA
            LB : c57bl6 fetal liver
            SM : c57bl6 mouse
            PL : ILLUMINA
            PM : HiSeq
            CN : CGSB
            Concentration : 0.475
            Barcode       : GAACTGAGCG-CGCTCCACGA
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Barcode @GACGAGATTAAGGATAATGT
            ID : CBJLFACXX:GACGAGATTAAGGATAATGT
            PU : CBJLFACXX:GACGAGATTAAGGATAATGT
            LB : c57bl6 adult bone marrow
            SM : c57bl6 mouse
            PL : ILLUMINA
            PM : HiSeq
            CN : CGSB
            Concentration : 0.475
            Barcode       : GACGAGATTA-AGGATAATGT
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Output feed No.0
            Type : sam
            Resolution : 2
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 4096
            URL : /dev/stdout?format=sam&compression=none



While not strictly necessary, You may examine the [compile configuration](example/CBJLFACXX_l01_sample_compiled.json), which is the actual configuration Pheniqs will execute will all implicit and default parameters. This is an easy way to see exactly what Pheniqs will be doing and spotting any configuration errors.

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json --compile
```

## Output format

By default Pheniqs will output uncompressed SAM to standard output. This allows you to quickly examine the output.

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json|less
```

You can see Pheniqs will declare the annotated read groups and populate the [BC](glossary.html#bc_auxiliary_tag), [QT](glossary.html#qt_auxiliary_tag) and [XB](glossary.html#xb_auxiliary_tag) tags.

    @HD     VN:1.0  SO:unknown      GO:query
    @RG     ID:CBJLFACXX:1:undetermined CN:CGSB PL:ILLUMINA PM:HiSeq  PU:CBJLFACXX:1:undetermined
    @RG     ID:CBJLFACXX:1:GAACTGAGCGCGCTCCACGA CN:CGSB PL:ILLUMINA PM:HiSeq PU:CBJLFACXX:1:GAACTGAGCGCGCTCCACGA SM:c57bl6 mouse LB:c57bl6 fetal liver
    @RG     ID:CBJLFACXX:1:GACGAGATTAAGGATAATGT CN:CGSB PL:ILLUMINA PM:HiSeq PU:CBJLFACXX:1:GACGAGATTAAGGATAATGT SM:c57bl6 mouse LB:c57bl6 adult bone marrow
    @PG     ID:pheniqs  PN:pheniqs  CL:pheniqs mux --config CBJLFACXX_l01_sample.json VN:2.0.6-70-g4093b3ab01c40c9093a4d95747a5a15518220e8a
    SN7001341:427:CBJLFACXX:1:1108:2455:1970 77  * 0 0 * * 0 0 TCCCTCTCACTGGAAGTGGTTGATCTCCAGGGAATCCCCAAGGTTAGCCT FFFF<B<00<BFFBB00<BF0BB0BFBBFF<<BBFB<FFFBBF707<BB RG:Z:CBJLFACXX:TCAGAC BC:Z:GACGACATTAAGGATAATGT QT:Z:IFFFIIIFFFFFBBFFFBFF XB:f:0.000293472
    SN7001341:427:CBJLFACXX:1:1108:2455:1970 141 * 0 0 * * 0 0 TTGATATTGTAAATTATAGACCAGACTGTGTACCATACTATTAATTGTCA <<'07'B<'<0'''0B'<<7<00<'<7'7<'0B7''<70<0<<7'<<'< RG:Z:CBJLFACXX:TCAGAC BC:Z:GACGACATTAAGGATAATGT QT:Z:IFFFIIIFFFFFBBFFFBFF XB:f:0.000293472

To write the from all multiplexed libraries into one BAM output file you can specify the `output` directive in the root of the configuration file
>```json
"output": [
    "CBJLFACXX_lane1_sample_demultiplex.bam"
],
```

Or directly on the command line

>```shell
pheniqs mux --config CBJLFACXX_l01_sample.json --output CBJLFACXX_lane1_sample_demultiplex.bam
```

If you want to split the libraries and their segments into separate fastq files you can specify the output on the individual barcode directives. Here are the [complete configuration](example/CBJLFACXX_l01_sample_split.json) and [compiled configuration](example/CBJLFACXX_l01_sample_split_compiled.json) for that scenario, but briefly those the changes:

>```json
"multiplex": {
  "codec": {
      "@GAACTGAGCGCGCTCCACGA": {
          "barcode": [ "GAACTGAGCG", "CGCTCCACGA" ],
          "output": [
              "CBJLFACXX_GAACTGAGCGCGCTCCACGA_l01s01.fastq.gz",
              "CBJLFACXX_GAACTGAGCGCGCTCCACGA_l01s02.fastq.gz"
          ]
      },
      "@GACGAGATTAAGGATAATGT": {
          "barcode": [ "GACGAGATTA", "AGGATAATGT" ],
          "output": [
              "CBJLFACXX_GACGAGATTAAGGATAATGT_l01s01.fastq.gz",
              "CBJLFACXX_GACGAGATTAAGGATAATGT_l01s02.fastq.gz"
          ]
      }
  },
  "undetermined": {
      "output": [
          "CBJLFACXX_undetermined_l01s01.fastq.gz",
          "CBJLFACXX_undetermined_l01s02.fastq.gz"
      ]
  }
}
```
>**Splitting the reads from different libraries** To write the segments of the different libraries to fastq file specify the output on the individual barcode directives. Notice we explicitly declare where undetermined reads will be written, otherwise they will be written to the default output, which is SAM on stdout.
{: .example}

## Prior estimation

Better estimation of the prior distribution of the samples can improve accuracy. Pheniqs provides a simple python script for adjusting your configuration to include priors estimated from the report emitted by a preliminary Pheniqs run. The `pheniqs-prior-api.py` script distributed with Pheniqs will execute Pheniqs with your given configuration and a special optimized mode that refrains from writing the output reads to save time and then emit a modified configuration file with adjusted priors. The priors you specify in your initial configuration can be your best guess for the priors but you can simply leave them out altogether.

>```shell
pheniqs-prior-api.py --configuration CBJLFACXX_l01_sample.json
```
