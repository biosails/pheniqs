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

This tutorial will walk you through hand writing configuration files for decoding sample barcodes in a standard Illumia sequencing run. For a tutorial that uses the python api to generate those configuration files see the [Standard Illumina sample decoding with the python API](illumina_python_api.html).

In this example the run has paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol. The read is made of 4 segments, 2 biological from the DNA or RNA fragment and 2 technical containing the i5 and i7 indices. If the results are written to SAM, BAM or CRAM the multiplex barcode and its quality scores are written to the [BC](glossary.html#bc_auxiliary_tag) and [QT](glossary.html#qt_auxiliary_tag) tags, and the decoding error probabilities is written to [XB](glossary.html#xb_auxiliary_tag).

## Input Read Layout

![paird end sequencing](/pheniqs/assets/img/paired_end_sequencing.png)

Base calling with bc2fastq will produce 4 files per lane: `H7LT2DSXX_S1_L001_R1_001.fastq.gz` containing the 3 prime prefix of the insert region, `H7LT2DSXX_S1_L001_I1_001.fastq.gz` containing the i7 index, `H7LT2DSXX_S1_L001_I2_001.fastq.gz` containing the i5 index and `H7LT2DSXX_S1_L001_R2_001.fastq.gz` containing the reverse complemented 5 prime suffix of the insert region, since it was read in reverse.

>```json
"input": [
    "H7LT2DSXX_S1_L001_R1_001.fastq.gz",
    "H7LT2DSXX_S1_L001_I1_001.fastq.gz",
    "H7LT2DSXX_S1_L001_I2_001.fastq.gz",
    "H7LT2DSXX_S1_L001_R2_001.fastq.gz"
],
```
>**declaring input read segments** 2 biological and 2 technical sequences are often found in 4 fastq files produced by bcl2fastq base calling.
{: .example}

To emit the two ends of the insert region as two segments of the output read we declare the global transform directives
>```json
"template": {
    "transform": { "token": [ "0::", "3::" ] }
}
```
>**declaring output read segments** Only the segments coming from the first and the forth file are biological sequences and should be included in the output.
{: .example}

To classify the reads by the i5 and i7 indices we declare the list of possible barcode sequences and a transform that tells Pheniqs where to find the barcode sequence.
>```json
"multiplex": {
    "transform": { "token": [ "1::8", "2::8" ] },
    "codec": {
        "@CGAGGCTGGTAAGGAG": {
            "barcode": [
                "CGAGGCTG",
                "GTAAGGAG"
            ]
        },
        "@AAGAGGCAACTGCATA": {
            "barcode": [
                "AAGAGGCA",
                "ACTGCATA"
            ]
        }
    }
}
```
>**declaring sample demultiplexing** The standard Illumina dual index protocol allows up to 96 unique dual 10 base barcodes in the i5 and i7 region, but for this example we only use 2 for brevity sake.
{: .example}

To discard reads that failed the internal Illumina sequencer chastity filter we instruct Pheniqs to filter incoming *qc fail* reads
>```json
"filter incoming qc fail": true
```

Putting things together we get a configuration for demultiplexing the reads in the three files with *PAMLD*. We expect **5%** of reads to be foreign DNA that should not classify to any of the barcodes, for instance spiked PhiX sequences, so we specify the noise prior with the `"noise": 0.05` directive. To mark as *qc fail* read where the posterior probability of correctly decoding the barcode is bellow **0.95** we specify `"confidence threshold": 0.95`.

>```json
{
    "CN": "CGSB",
    "PL": "ILLUMINA",
    "PM": "A00534",
    "filter incoming qc fail": true,
    "flowcell id": "H7LT2DSXX",
    "flowcell lane number": 1,
    "input": [
        "H7LT2DSXX_S1_L001_R1_001.fastq.gz",
        "H7LT2DSXX_S1_L001_I1_001.fastq.gz",
        "H7LT2DSXX_S1_L001_I2_001.fastq.gz",
        "H7LT2DSXX_S1_L001_R2_001.fastq.gz"
    ],
    "multiplex": {
        "algorithm": "pamld",
        "codec": {
            "@CGAGGCTGGTAAGGAG": {
                "barcode": [
                    "CGAGGCTG",
                    "GTAAGGAG"
                ]
            },
            "@AAGAGGCAACTGCATA": {
                "barcode": [
                    "AAGAGGCA",
                    "ACTGCATA"
                ]
            }
        },
        "confidence threshold": 0.95,
        "noise": 0.05,
        "transform": {
            "token": [
                "1::8",
                "2::8"
            ]
        }
    },
    "template": {
        "transform": {
            "token": [
                "0::",
                "3::"
            ]
        }
    }
}
```
>**Single index Paired end Illumina protocol** Classifying the 2 barcodes using the 4 fastq files produced by bcl2fastq.
{: .example}

Before we proceed we validate the configuration with Pheniqs

>```shell
pheniqs mux --config H7LT2DSXX_l01_sample.json --validate
```

The output is a readable description of all the explicit and implicit parameters after applying defaults. You can check how Pheniqs detects the input format, compression and layout as well as the output you can expect. In The *Output transform* section is a verbal description of how the output read segments will be assembled from the input. similarly *Transform* in the *Mutliplex decoding* section desribes how the segment that will be matched against the barcodes is assembled. the You can also see how each of the read groups will be tagged and the prior probability PAMLD will assume for each barcode.

    Environment

        Base input URL                              /home/lg/H7LT2DSXX
        Base output URL                             /home/lg/H7LT2DSXX
        Platform                                    ILLUMINA
        Quality tracking                            disabled
        Filter incoming QC failed reads             enabled
        Filter outgoing QC failed reads             disabled
        Input Phred offset                          33
        Output Phred offset                         33
        Leading segment index                       0
        Default output format                       sam
        Default output compression                  unknown
        Default output compression level            5
        Feed buffer capacity                        2048
        Threads                                     8
        Decoding threads                            1
        HTSLib threads                              8

    Input

        Input segment cardinality                   4

        Input segment No.0 : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz
        Input segment No.1 : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz
        Input segment No.2 : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz
        Input segment No.3 : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz

        Input feed No.0
            Type : fastq
            Compression : gz
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz

        Input feed No.1
            Type : fastq
            Compression : gz
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz

        Input feed No.2
            Type : fastq
            Compression : gz
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz

        Input feed No.3
            Type : fastq
            Compression : gz
            Resolution : 1
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 2048
            URL : /home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz

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
        Shannon bound                               3 3
        Noise                                       0.05
        Confidence threshold                        0.95
        Segment cardinality                         2
        Nucleotide cardinality                      16
        Barcode segment length                      8 8

        Transform

            Token No.0
                Length        8
                Pattern       1::8
                Description   cycles 0 to 8 of input segment 1

            Token No.1
                Length        8
                Pattern       2::8
                Description   cycles 0 to 8 of input segment 2

            Assembly instruction
                Append token 0 of input segment 1 to output segment 0
                Append token 1 of input segment 2 to output segment 1


        Barcode undetermined
            ID : H7LT2DSXX:1:undetermined
            PU : H7LT2DSXX:1:undetermined
            PL : ILLUMINA
            PM : A00534
            CN : CGSB
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Barcode @AAGAGGCAACTGCATA
            ID : H7LT2DSXX:1:AAGAGGCAACTGCATA
            PU : H7LT2DSXX:1:AAGAGGCAACTGCATA
            PL : ILLUMINA
            PM : A00534
            CN : CGSB
            Concentration : 0.475
            Barcode       : AAGAGGCA-ACTGCATA
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Barcode @CGAGGCTGGTAAGGAG
            ID : H7LT2DSXX:1:CGAGGCTGGTAAGGAG
            PU : H7LT2DSXX:1:CGAGGCTGGTAAGGAG
            PL : ILLUMINA
            PM : A00534
            CN : CGSB
            Concentration : 0.475
            Barcode       : CGAGGCTG-GTAAGGAG
            Segment No.0  : /dev/stdout?format=sam&compression=none
            Segment No.1  : /dev/stdout?format=sam&compression=none

        Output feed No.0
            Type : sam
            Resolution : 2
            Phred offset : 33
            Platform : ILLUMINA
            Buffer capacity : 4096
            URL : /dev/stdout?format=sam&compression=none

While not strictly necessary, You may examine the [compile configuration](example/H7LT2DSXX_l01_sample_compiled.json), which is the actual configuration Pheniqs will execute will all implicit and default parameters. This is an easy way to see exactly what Pheniqs will be doing and spotting any configuration errors.

>```shell
pheniqs mux --config H7LT2DSXX_l01_sample.json --compile
```

>```json
{
    "CN": "CGSB",
    "PL": "ILLUMINA",
    "PM": "A00534",
    "application name": "pheniqs",
    "base input url": "/home/lg/H7LT2DSXX",
    "base output url": "/home/lg/H7LT2DSXX",
    "buffer capacity": 2048,
    "decoding threads": 1,
    "default output compression": "unknown",
    "default output compression level": "5",
    "default output format": "sam",
    "feed": {
        "input feed": [
            {
                "capacity": 2048,
                "direction": "in",
                "index": 0,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 1,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 2,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 3,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz"
            }
        ],
        "input feed by segment": [
            {
                "capacity": 2048,
                "direction": "in",
                "index": 0,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 1,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 2,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz"
            },
            {
                "capacity": 2048,
                "direction": "in",
                "index": 3,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 1,
                "url": "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz"
            }
        ],
        "output feed": [
            {
                "capacity": 4096,
                "direction": "out",
                "index": 0,
                "phred offset": 33,
                "platform": "ILLUMINA",
                "resolution": 2,
                "url": "/dev/stdout?format=sam&compression=none"
            }
        ]
    },
    "filter incoming qc fail": true,
    "float precision": 15,
    "flowcell id": "H7LT2DSXX",
    "flowcell lane number": 1,
    "full command": "pheniqs mux --config /Users/lg/Desktop/tmp --compile",
    "htslib threads": 8,
    "input": [
        "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R1_001.fastq.gz?format=fastq&compression=gz",
        "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I1_001.fastq.gz?format=fastq&compression=gz",
        "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_I2_001.fastq.gz?format=fastq&compression=gz",
        "/home/lg/H7LT2DSXX/H7LT2DSXX_S1_L001_R2_001.fastq.gz?format=fastq&compression=gz"
    ],
    "input phred offset": 33,
    "input segment cardinality": 4,
    "leading segment index": 0,
    "multiplex": {
        "CN": "CGSB",
        "PL": "ILLUMINA",
        "PM": "A00534",
        "algorithm": "pamld",
        "barcode cardinality": 3,
        "barcode length": [
            8,
            8
        ],
        "base output url": "/home/lg/H7LT2DSXX",
        "codec": {
            "@AAGAGGCAACTGCATA": {
                "CN": "CGSB",
                "ID": "H7LT2DSXX:1:AAGAGGCAACTGCATA",
                "PL": "ILLUMINA",
                "PM": "A00534",
                "PU": "H7LT2DSXX:1:AAGAGGCAACTGCATA",
                "TC": 2,
                "algorithm": "pamld",
                "barcode": [
                    "AAGAGGCA",
                    "ACTGCATA"
                ],
                "concentration": 0.475,
                "feed": {
                    "output feed by segment": [
                        {
                            "capacity": 4096,
                            "direction": "out",
                            "index": 0,
                            "phred offset": 33,
                            "platform": "ILLUMINA",
                            "resolution": 2,
                            "url": "/dev/stdout?format=sam&compression=none"
                        },
                        {
                            "capacity": 4096,
                            "direction": "out",
                            "index": 0,
                            "phred offset": 33,
                            "platform": "ILLUMINA",
                            "resolution": 2,
                            "url": "/dev/stdout?format=sam&compression=none"
                        }
                    ]
                },
                "flowcell id": "H7LT2DSXX",
                "flowcell lane number": 1,
                "index": 1,
                "output": [
                    "/dev/stdout?format=sam&compression=none",
                    "/dev/stdout?format=sam&compression=none"
                ],
                "segment cardinality": 2
            },
            "@CGAGGCTGGTAAGGAG": {
                "CN": "CGSB",
                "ID": "H7LT2DSXX:1:CGAGGCTGGTAAGGAG",
                "PL": "ILLUMINA",
                "PM": "A00534",
                "PU": "H7LT2DSXX:1:CGAGGCTGGTAAGGAG",
                "TC": 2,
                "algorithm": "pamld",
                "barcode": [
                    "CGAGGCTG",
                    "GTAAGGAG"
                ],
                "concentration": 0.475,
                "feed": {
                    "output feed by segment": [
                        {
                            "capacity": 4096,
                            "direction": "out",
                            "index": 0,
                            "phred offset": 33,
                            "platform": "ILLUMINA",
                            "resolution": 2,
                            "url": "/dev/stdout?format=sam&compression=none"
                        },
                        {
                            "capacity": 4096,
                            "direction": "out",
                            "index": 0,
                            "phred offset": 33,
                            "platform": "ILLUMINA",
                            "resolution": 2,
                            "url": "/dev/stdout?format=sam&compression=none"
                        }
                    ]
                },
                "flowcell id": "H7LT2DSXX",
                "flowcell lane number": 1,
                "index": 2,
                "output": [
                    "/dev/stdout?format=sam&compression=none",
                    "/dev/stdout?format=sam&compression=none"
                ],
                "segment cardinality": 2
            }
        },
        "confidence threshold": 0.95,
        "distance tolerance": [
            3,
            3
        ],
        "flowcell id": "H7LT2DSXX",
        "flowcell lane number": 1,
        "index": 0,
        "noise": 0.05,
        "nucleotide cardinality": 16,
        "output": [
            "/dev/stdout"
        ],
        "quality masking threshold": 0,
        "segment cardinality": 2,
        "shannon bound": [
            3,
            3
        ],
        "transform": {
            "knit": [
                "0",
                "1"
            ],
            "token": [
                "1::8",
                "2::8"
            ]
        },
        "undetermined": {
            "CN": "CGSB",
            "ID": "H7LT2DSXX:1:undetermined",
            "PL": "ILLUMINA",
            "PM": "A00534",
            "PU": "H7LT2DSXX:1:undetermined",
            "TC": 2,
            "algorithm": "pamld",
            "barcode": [
                "========",
                "========"
            ],
            "concentration": 0.05,
            "feed": {
                "output feed by segment": [
                    {
                        "capacity": 4096,
                        "direction": "out",
                        "index": 0,
                        "phred offset": 33,
                        "platform": "ILLUMINA",
                        "resolution": 2,
                        "url": "/dev/stdout?format=sam&compression=none"
                    },
                    {
                        "capacity": 4096,
                        "direction": "out",
                        "index": 0,
                        "phred offset": 33,
                        "platform": "ILLUMINA",
                        "resolution": 2,
                        "url": "/dev/stdout?format=sam&compression=none"
                    }
                ]
            },
            "flowcell id": "H7LT2DSXX",
            "flowcell lane number": 1,
            "index": 0,
            "output": [
                "/dev/stdout?format=sam&compression=none",
                "/dev/stdout?format=sam&compression=none"
            ],
            "segment cardinality": 2
        }
    },
    "output": [
        "/dev/stdout"
    ],
    "output phred offset": 33,
    "output segment cardinality": 2,
    "platform": "ILLUMINA",
    "report url": "/dev/stderr",
    "template": {
        "transform": {
            "knit": [
                "0",
                "1"
            ],
            "token": [
                "0::",
                "3::"
            ]
        }
    },
    "threads": 8
}
```
>**Compiled configuration** compiling a configuration is a good way to see what all the implicit parameter values are going to be.
{: .example}


## Output format

By default Pheniqs will output uncompressed SAM to standard output. This allows you to quickly examine the output. To write the from all multiplexed libraries into one BAM output file you can specify the `output` directive in the root of the configuration file

>```json
"output": [
    "H7LT2DSXX_l01_sample.bam"
],
```
>**adding an output directive** by default output is written to stdout.
{: .example}

Or directly on the command line

>```shell
pheniqs mux --config H7LT2DSXX_l01_sample.json --output H7LT2DSXX_l01_sample.bam
```
