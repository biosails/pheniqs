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
        <li><a                  href="/pheniqs/">Home</a></li>
        <li><a                  href="/pheniqs/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a class="active"   href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Workflow
{:.page-title}

## Complete examples

```json
{
    "CN": "CGSB", 
    "DT": "2016-07-15T07:00:00+00:00", 
    "PI": "1500", 
    "PL": "ILLUMINA", 
    "PM": "hiseq 2500", 
    "base input path": "~/HK5NHBGXX/1", 
    "base output path": "~/HK5NHBGXX/1", 
    "channel": [
        {
            "DS": "undetermined description", 
            "LB": "undetermined_library", 
            "PU": "HK5NHBGXX:1:undetermined", 
            "RG": "HK5NHBGXX:1:undetermined", 
            "SM": "undetermined_sample", 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ], 
            "undetermined": true
        }, 
        {
            "DS": "ntr_2 description", 
            "LB": "ntr_2_library", 
            "PU": "HK5NHBGXX:1:CGTACTAGTCTTACGC", 
            "RG": "HK5NHBGXX:1:CGTACTAGTCTTACGC", 
            "SM": "ntr_2_sample", 
            "barcode": [
                "CGTACTAG", 
                "TCTTACGC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_10 description", 
            "LB": "ntr_10_library", 
            "PU": "HK5NHBGXX:1:CGAGGCTGTCTTACGC", 
            "RG": "HK5NHBGXX:1:CGAGGCTGTCTTACGC", 
            "SM": "ntr_10_sample", 
            "barcode": [
                "CGAGGCTG", 
                "TCTTACGC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_11 description", 
            "LB": "ntr_11_library", 
            "PU": "HK5NHBGXX:1:AAGAGGCATCTTACGC", 
            "RG": "HK5NHBGXX:1:AAGAGGCATCTTACGC", 
            "SM": "ntr_11_sample", 
            "barcode": [
                "AAGAGGCA", 
                "TCTTACGC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_12 description", 
            "LB": "ntr_12_library", 
            "PU": "HK5NHBGXX:1:GTAGAGGATCTTACGC", 
            "RG": "HK5NHBGXX:1:GTAGAGGATCTTACGC", 
            "SM": "ntr_12_sample", 
            "barcode": [
                "GTAGAGGA", 
                "TCTTACGC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_13 description", 
            "LB": "ntr_13_library", 
            "PU": "HK5NHBGXX:1:TAAGGCGAATAGAGAG", 
            "RG": "HK5NHBGXX:1:TAAGGCGAATAGAGAG", 
            "SM": "ntr_13_sample", 
            "barcode": [
                "TAAGGCGA", 
                "ATAGAGAG"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_14 description", 
            "LB": "ntr_14_library", 
            "PU": "HK5NHBGXX:1:CGTACTAGATAGAGAG", 
            "RG": "HK5NHBGXX:1:CGTACTAGATAGAGAG", 
            "SM": "ntr_14_sample", 
            "barcode": [
                "CGTACTAG", 
                "ATAGAGAG"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_36 description", 
            "LB": "ntr_36_library", 
            "PU": "HK5NHBGXX:1:GTAGAGGAAGAGGATA", 
            "RG": "HK5NHBGXX:1:GTAGAGGAAGAGGATA", 
            "SM": "ntr_36_sample", 
            "barcode": [
                "GTAGAGGA", 
                "AGAGGATA"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_37 description", 
            "LB": "ntr_37_library", 
            "PU": "HK5NHBGXX:1:TAAGGCGATCTACTCT", 
            "RG": "HK5NHBGXX:1:TAAGGCGATCTACTCT", 
            "SM": "ntr_37_sample", 
            "barcode": [
                "TAAGGCGA", 
                "TCTACTCT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_58 description", 
            "LB": "ntr_58_library", 
            "PU": "HK5NHBGXX:1:CGAGGCTGCTCCTTAC", 
            "RG": "HK5NHBGXX:1:CGAGGCTGCTCCTTAC", 
            "SM": "ntr_58_sample", 
            "barcode": [
                "CGAGGCTG", 
                "CTCCTTAC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_59 description", 
            "LB": "ntr_59_library", 
            "PU": "HK5NHBGXX:1:AAGAGGCACTCCTTAC", 
            "RG": "HK5NHBGXX:1:AAGAGGCACTCCTTAC", 
            "SM": "ntr_59_sample", 
            "barcode": [
                "AAGAGGCA", 
                "CTCCTTAC"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_67 description", 
            "LB": "ntr_67_library", 
            "PU": "HK5NHBGXX:1:CTCTCTACTATGCAGT", 
            "RG": "HK5NHBGXX:1:CTCTCTACTATGCAGT", 
            "SM": "ntr_67_sample", 
            "barcode": [
                "CTCTCTAC", 
                "TATGCAGT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_68 description", 
            "LB": "ntr_68_library", 
            "PU": "HK5NHBGXX:1:CAGAGAGGTATGCAGT", 
            "RG": "HK5NHBGXX:1:CAGAGAGGTATGCAGT", 
            "SM": "ntr_68_sample", 
            "barcode": [
                "CAGAGAGG", 
                "TATGCAGT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_69 description", 
            "LB": "ntr_69_library", 
            "PU": "HK5NHBGXX:1:GCTACGCTTATGCAGT", 
            "RG": "HK5NHBGXX:1:GCTACGCTTATGCAGT", 
            "SM": "ntr_69_sample", 
            "barcode": [
                "GCTACGCT", 
                "TATGCAGT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_76 description", 
            "LB": "ntr_76_library", 
            "PU": "HK5NHBGXX:1:TCCTGAGCTACTCCTT", 
            "RG": "HK5NHBGXX:1:TCCTGAGCTACTCCTT", 
            "SM": "ntr_76_sample", 
            "barcode": [
                "TCCTGAGC", 
                "TACTCCTT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_77 description", 
            "LB": "ntr_77_library", 
            "PU": "HK5NHBGXX:1:GGACTCCTTACTCCTT", 
            "RG": "HK5NHBGXX:1:GGACTCCTTACTCCTT", 
            "SM": "ntr_77_sample", 
            "barcode": [
                "GGACTCCT", 
                "TACTCCTT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_78 description", 
            "LB": "ntr_78_library", 
            "PU": "HK5NHBGXX:1:TAGGCATGTACTCCTT", 
            "RG": "HK5NHBGXX:1:TAGGCATGTACTCCTT", 
            "SM": "ntr_78_sample", 
            "barcode": [
                "TAGGCATG", 
                "TACTCCTT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_83 description", 
            "LB": "ntr_83_library", 
            "PU": "HK5NHBGXX:1:AAGAGGCATACTCCTT", 
            "RG": "HK5NHBGXX:1:AAGAGGCATACTCCTT", 
            "SM": "ntr_83_sample", 
            "barcode": [
                "AAGAGGCA", 
                "TACTCCTT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_84 description", 
            "LB": "ntr_84_library", 
            "PU": "HK5NHBGXX:1:GTAGAGGATACTCCTT", 
            "RG": "HK5NHBGXX:1:GTAGAGGATACTCCTT", 
            "SM": "ntr_84_sample", 
            "barcode": [
                "GTAGAGGA", 
                "TACTCCTT"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_87 description", 
            "LB": "ntr_87_library", 
            "PU": "HK5NHBGXX:1:AGGCAGAAAGGCTTAG", 
            "RG": "HK5NHBGXX:1:AGGCAGAAAGGCTTAG", 
            "SM": "ntr_87_sample", 
            "barcode": [
                "AGGCAGAA", 
                "AGGCTTAG"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_88 description", 
            "LB": "ntr_88_library", 
            "PU": "HK5NHBGXX:1:TCCTGAGCAGGCTTAG", 
            "RG": "HK5NHBGXX:1:TCCTGAGCAGGCTTAG", 
            "SM": "ntr_88_sample", 
            "barcode": [
                "TCCTGAGC", 
                "AGGCTTAG"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }, 
        {
            "DS": "ntr_96 description", 
            "LB": "ntr_96_library", 
            "PU": "HK5NHBGXX:1:GTAGAGGAAGGCTTAG", 
            "RG": "HK5NHBGXX:1:GTAGAGGAAGGCTTAG", 
            "SM": "ntr_96_sample", 
            "barcode": [
                "GTAGAGGA", 
                "AGGCTTAG"
            ], 
            "concentration": 1, 
            "output": [
                "HK5NHBGXX_l01n01.cram"
            ]
        }
    ], 
    "confidence": 0.99, 
    "decoder": "pamld", 
    "include filtered": false, 
    "input": [
        "HK5NHBGXX_1_interleaved.fastq.gz", 
        "HK5NHBGXX_1_interleaved.fastq.gz", 
        "HK5NHBGXX_1_interleaved.fastq.gz", 
        "HK5NHBGXX_1_interleaved.fastq.gz"
    ], 
    "molecular barcode": [
        "~4"
    ], 
    "multiplex barcode": [
        "~2", 
        "3"
    ], 
    "noise": 0.01, 
    "threads": 4, 
    "template": [
        "0", 
        "1"
    ], 
    "token": [
        "0:0:", 
        "3:0:", 
        "1:0:8", 
        "2:0:8", 
        "0:-7:-1"
    ]
}
```

To validate a configuration you can execute

```zsh
pheniqs demux --config HK5NHBGXX.json --validate
```

Which prints the following validation report to standard output

```
Environment 

    Version                                     0.9.d55d65203560ae3719fcbb79ac657189e463ffed
    Barcode decoder                             pamld
    Platform                                    ILLUMINA
    Long read                                   disabled
    Quality tracking                            enabled
    Leading template segment                    0
    Input Phred offset                          33
    Output Phred offset                         33
    Include non PF reads                        disabled
    Undetermined reads                          Channel No.0
    Multiplex barcoding cycles                  16
    Shortest multiplex distance                 5
    Multiplex barcode length                    8 8 
    Prior noise frequency                       0.010000
    Multiplex confidence                        0.990000
    Random word frequency                       0.0000000002328306437
    Molecular barcode length                    6 
    Feed buffer capacity                        2048
    Transforming threads                        1
    Threads                                     1

Transformation 

    Token No.0
        Length        variable
        Pattern       0::
        Description   cycles 0 to end of input segment 0

    Token No.1
        Length        variable
        Pattern       3::
        Description   cycles 0 to end of input segment 3

    Token No.2
        Length        8
        Pattern       1::8
        Description   cycles 0 to 8 of input segment 1

    Token No.3
        Length        8
        Pattern       2::8
        Description   cycles 0 to 8 of input segment 2

    Token No.4
        Length        6
        Pattern       0:-7:-1
        Description   cycles -7 to -1 of input segment 0

    Template transform
        Append token 0 of input segment 0 to output segment 0
        Append token 1 of input segment 3 to output segment 1

    Multiplex barcode transform
        Append reverse complemented token 2 of input segment 1 to output segment 0
        Append token 3 of input segment 2 to output segment 1

    Molecular barcode transform
        Append reverse complemented token 4 of input segment 0 to output segment 0

Input 

    Input segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_1_interleaved.fastq.gz
    Input segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_1_interleaved.fastq.gz
    Input segment No.2 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_1_interleaved.fastq.gz
    Input segment No.3 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_1_interleaved.fastq.gz

Channel 

    Channel No.0
        RG : HK5NHBGXX:1:undetermined
        PU : HK5NHBGXX:1:undetermined
        LB : undetermined_library
        SM : undetermined_sample
        CN : CGSB
        DS : undetermined description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        Undetermined : true

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.1
        RG : HK5NHBGXX:1:CGTACTAGTCTTACGC
        PU : HK5NHBGXX:1:CGTACTAGTCTTACGC
        LB : ntr_2_library
        SM : ntr_2_sample
        CN : CGSB
        DS : ntr_2 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CGTACTAG
        Multiplex barcode No.1 : TCTTACGC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.2
        RG : HK5NHBGXX:1:CGAGGCTGTCTTACGC
        PU : HK5NHBGXX:1:CGAGGCTGTCTTACGC
        LB : ntr_10_library
        SM : ntr_10_sample
        CN : CGSB
        DS : ntr_10 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CGAGGCTG
        Multiplex barcode No.1 : TCTTACGC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.3
        RG : HK5NHBGXX:1:AAGAGGCATCTTACGC
        PU : HK5NHBGXX:1:AAGAGGCATCTTACGC
        LB : ntr_11_library
        SM : ntr_11_sample
        CN : CGSB
        DS : ntr_11 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : AAGAGGCA
        Multiplex barcode No.1 : TCTTACGC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.4
        RG : HK5NHBGXX:1:GTAGAGGATCTTACGC
        PU : HK5NHBGXX:1:GTAGAGGATCTTACGC
        LB : ntr_12_library
        SM : ntr_12_sample
        CN : CGSB
        DS : ntr_12 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GTAGAGGA
        Multiplex barcode No.1 : TCTTACGC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.5
        RG : HK5NHBGXX:1:TAAGGCGAATAGAGAG
        PU : HK5NHBGXX:1:TAAGGCGAATAGAGAG
        LB : ntr_13_library
        SM : ntr_13_sample
        CN : CGSB
        DS : ntr_13 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : TAAGGCGA
        Multiplex barcode No.1 : ATAGAGAG

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.6
        RG : HK5NHBGXX:1:CGTACTAGATAGAGAG
        PU : HK5NHBGXX:1:CGTACTAGATAGAGAG
        LB : ntr_14_library
        SM : ntr_14_sample
        CN : CGSB
        DS : ntr_14 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CGTACTAG
        Multiplex barcode No.1 : ATAGAGAG

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.7
        RG : HK5NHBGXX:1:GTAGAGGAAGAGGATA
        PU : HK5NHBGXX:1:GTAGAGGAAGAGGATA
        LB : ntr_36_library
        SM : ntr_36_sample
        CN : CGSB
        DS : ntr_36 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GTAGAGGA
        Multiplex barcode No.1 : AGAGGATA

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.8
        RG : HK5NHBGXX:1:TAAGGCGATCTACTCT
        PU : HK5NHBGXX:1:TAAGGCGATCTACTCT
        LB : ntr_37_library
        SM : ntr_37_sample
        CN : CGSB
        DS : ntr_37 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : TAAGGCGA
        Multiplex barcode No.1 : TCTACTCT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.9
        RG : HK5NHBGXX:1:CGAGGCTGCTCCTTAC
        PU : HK5NHBGXX:1:CGAGGCTGCTCCTTAC
        LB : ntr_58_library
        SM : ntr_58_sample
        CN : CGSB
        DS : ntr_58 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CGAGGCTG
        Multiplex barcode No.1 : CTCCTTAC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.10
        RG : HK5NHBGXX:1:AAGAGGCACTCCTTAC
        PU : HK5NHBGXX:1:AAGAGGCACTCCTTAC
        LB : ntr_59_library
        SM : ntr_59_sample
        CN : CGSB
        DS : ntr_59 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : AAGAGGCA
        Multiplex barcode No.1 : CTCCTTAC

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.11
        RG : HK5NHBGXX:1:CTCTCTACTATGCAGT
        PU : HK5NHBGXX:1:CTCTCTACTATGCAGT
        LB : ntr_67_library
        SM : ntr_67_sample
        CN : CGSB
        DS : ntr_67 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CTCTCTAC
        Multiplex barcode No.1 : TATGCAGT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.12
        RG : HK5NHBGXX:1:CAGAGAGGTATGCAGT
        PU : HK5NHBGXX:1:CAGAGAGGTATGCAGT
        LB : ntr_68_library
        SM : ntr_68_sample
        CN : CGSB
        DS : ntr_68 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : CAGAGAGG
        Multiplex barcode No.1 : TATGCAGT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.13
        RG : HK5NHBGXX:1:GCTACGCTTATGCAGT
        PU : HK5NHBGXX:1:GCTACGCTTATGCAGT
        LB : ntr_69_library
        SM : ntr_69_sample
        CN : CGSB
        DS : ntr_69 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GCTACGCT
        Multiplex barcode No.1 : TATGCAGT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.14
        RG : HK5NHBGXX:1:TCCTGAGCTACTCCTT
        PU : HK5NHBGXX:1:TCCTGAGCTACTCCTT
        LB : ntr_76_library
        SM : ntr_76_sample
        CN : CGSB
        DS : ntr_76 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : TCCTGAGC
        Multiplex barcode No.1 : TACTCCTT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.15
        RG : HK5NHBGXX:1:GGACTCCTTACTCCTT
        PU : HK5NHBGXX:1:GGACTCCTTACTCCTT
        LB : ntr_77_library
        SM : ntr_77_sample
        CN : CGSB
        DS : ntr_77 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GGACTCCT
        Multiplex barcode No.1 : TACTCCTT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.16
        RG : HK5NHBGXX:1:TAGGCATGTACTCCTT
        PU : HK5NHBGXX:1:TAGGCATGTACTCCTT
        LB : ntr_78_library
        SM : ntr_78_sample
        CN : CGSB
        DS : ntr_78 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : TAGGCATG
        Multiplex barcode No.1 : TACTCCTT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.17
        RG : HK5NHBGXX:1:AAGAGGCATACTCCTT
        PU : HK5NHBGXX:1:AAGAGGCATACTCCTT
        LB : ntr_83_library
        SM : ntr_83_sample
        CN : CGSB
        DS : ntr_83 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : AAGAGGCA
        Multiplex barcode No.1 : TACTCCTT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.18
        RG : HK5NHBGXX:1:GTAGAGGATACTCCTT
        PU : HK5NHBGXX:1:GTAGAGGATACTCCTT
        LB : ntr_84_library
        SM : ntr_84_sample
        CN : CGSB
        DS : ntr_84 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GTAGAGGA
        Multiplex barcode No.1 : TACTCCTT

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.19
        RG : HK5NHBGXX:1:AGGCAGAAAGGCTTAG
        PU : HK5NHBGXX:1:AGGCAGAAAGGCTTAG
        LB : ntr_87_library
        SM : ntr_87_sample
        CN : CGSB
        DS : ntr_87 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : AGGCAGAA
        Multiplex barcode No.1 : AGGCTTAG

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.20
        RG : HK5NHBGXX:1:TCCTGAGCAGGCTTAG
        PU : HK5NHBGXX:1:TCCTGAGCAGGCTTAG
        LB : ntr_88_library
        SM : ntr_88_sample
        CN : CGSB
        DS : ntr_88 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : TCCTGAGC
        Multiplex barcode No.1 : AGGCTTAG

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

    Channel No.21
        RG : HK5NHBGXX:1:GTAGAGGAAGGCTTAG
        PU : HK5NHBGXX:1:GTAGAGGAAGGCTTAG
        LB : ntr_96_library
        SM : ntr_96_sample
        CN : CGSB
        DS : ntr_96 description
        DT : 2016-07-15T07:00:00+00:00
        PL : ILLUMINA
        PM : hiseq 2500
        PI : 1500
        TC : 2
        PC : 0.0471428571428571

        Multiplex barcode No.0 : GTAGAGGA
        Multiplex barcode No.1 : AGGCTTAG

        Output segment No.0 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram
        Output segment No.1 : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

Feed 

    Input feed No.0
        Type : fastq
        Compression : gz
        Resolution : 4
        Phred offset : 33
        Platform : ILLUMINA
        Buffer capacity : 8192
        URL : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_1_interleaved.fastq.gz

    Output feed No.0
        Type : cram
        Resolution : 2
        Phred offset : 33
        Platform : ILLUMINA
        Buffer capacity : 4096
        URL : /Users/lg/HK5NHBGXX/1/HK5NHBGXX_l01n01.cram

Barcode distance distribution

    0 4 4 7 7 8 5 7 6 5 6 7 AAGAGGCA 2 
    1 0 7 7 6 7 7 6 5 7 7 8 AGGCAGAA 2 
    1 3 0 5 5 7 7 8 7 5 4 5 CAGAGAGG 2 
    3 3 2 0 5 7 8 5 6 4 5 7 CGAGGCTG 2 
    3 2 2 2 0 4 5 7 8 8 6 8 CGTACTAG 2 
    3 3 3 3 1 0 7 8 7 8 7 5 CTCTCTAC 3 
    2 3 3 3 2 3 0 5 6 8 7 7 GCTACGCT 3 
    3 2 3 2 3 3 2 0 6 6 8 8 GGACTCCT 3 
    2 2 3 2 3 3 2 2 0 4 7 7 GTAGAGGA 2 
    2 3 2 1 3 3 3 2 1 0 5 5 TAAGGCGA 2 
    2 3 1 2 2 3 3 3 3 2 0 6 TAGGCATG 2 
    3 3 2 3 3 2 3 3 3 2 2 0 TCCTGAGC 3 

    0 6 5 8 7 6 8 8 AGAGGATA 3 
    2 0 5 4 8 8 7 8 AGGCTTAG 2 
    2 2 0 6 8 7 8 7 ATAGAGAG 2 
    3 1 2 0 7 8 7 7 CTCCTTAC 2 
    3 3 3 3 0 4 5 5 TACTCCTT 2 
    2 3 3 3 1 0 4 5 TATGCAGT 2 
    3 3 3 3 2 1 0 5 TCTACTCT 2 
    3 3 3 3 2 2 2 0 TCTTACGC 2 

    0  7  7  8  12 7  14 13 14 16 13 14 14 10 13 13 11 12 13 11 14 AAGAGGCACTCCTTAC 5  
    3  0  5  12 8  14 12 15 12 12 9  7  13 14 6  11 13 10 6  15 7  AAGAGGCATACTCCTT 4  
    3  2  0  12 9  14 7  14 7  13 10 12 14 14 11 6  12 10 11 15 12 AAGAGGCATCTTACGC 5  
    3  5  5  0  15 11 15 11 14 15 15 14 11 5  13 13 12 14 15 8  16 AGGCAGAAAGGCTTAG 5  
    5  3  4  7  0  13 10 12 10 7  7  12 13 15 11 12 12 9  8  13 9  CAGAGAGGTATGCAGT 5  
    3  6  6  5  6  0  7  11 12 15 16 12 14 10 13 13 10 11 12 11 14 CGAGGCTGCTCCTTAC 5  
    6  5  3  7  4  3  0  12 5  12 13 10 14 14 11 6  11 9  10 15 12 CGAGGCTGTCTTACGC 5  
    6  7  6  5  5  5  5  0  7  11 12 15 13 13 16 15 8  16 14 13 16 CGTACTAGATAGAGAG 6  
    6  5  3  6  4  5  2  3  0  9  10 12 16 16 13 8  15 13 11 16 13 CGTACTAGTCTTACGC 5  
    7  5  6  7  3  7  5  5  4  0  7  12 13 15 11 12 15 12 11 13 9  CTCTCTACTATGCAGT 5  
    6  4  4  7  3  7  6  5  4  3  0  9  12 14 10 11 15 12 11 15 11 GCTACGCTTATGCAGT 5  
    6  3  5  6  5  5  4  7  5  5  4  0  13 14 6  11 14 11 8  16 8  GGACTCCTTACTCCTT 5  
    6  6  6  5  6  6  6  6  7  6  5  6  0  6  7  8  9  12 14 13 14 GTAGAGGAAGAGGATA 5  
    4  6  6  2  7  4  6  6  7  7  6  6  2  0  8  8  9  11 15 7  15 GTAGAGGAAGGCTTAG 5  
    6  2  5  6  5  6  5  7  6  5  4  2  3  3  0  5  12 9  7  15 7  GTAGAGGATACTCCTT 4  
    6  5  2  6  5  6  2  7  3  5  5  5  3  3  2  0  11 9  12 15 12 GTAGAGGATCTTACGC 5  
    5  6  5  5  5  4  5  3  7  7  7  6  4  4  5  5  0  8  13 10 13 TAAGGCGAATAGAGAG 5  
    5  4  4  6  4  5  4  7  6  5  5  5  5  5  4  4  3  0  10 12 10 TAAGGCGATCTACTCT 5  
    6  2  5  7  3  5  4  6  5  5  5  3  6  7  3  5  6  4  0  14 6  TAGGCATGTACTCCTT 5  
    5  7  7  3  6  5  7  6  7  6  7  7  6  3  7  7  4  5  6  0  8  TCCTGAGCAGGCTTAG 6  
    6  3  5  7  4  6  5  7  6  4  5  3  6  7  3  5  6  4  2  3  0  TCCTGAGCTACTCCTT 5  
```