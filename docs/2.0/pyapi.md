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
        <li><a                  href="/pheniqs/2.0/pyapi.html">Python API</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Pheniqs python API
{:.page-title}

## Illumina API

`pheniqs-illumina-api.py` can generate Pheniqs configuration files from metadata found in an Illumina run folder.


>```shell
    usage: pheniqs-illumina-api.py [-h] [--version] [-v LEVEL] ACTION ...

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            logging verbosity level

    pipeline operations:
      Generate pheniqs configuration files or a bcl2fastq command from an
      illumina run directory. This tool parses that RunInfo.xml,
      RunParameters.xml and SampleSheet.csv files in the directory.

      ACTION
        basecall            bcl2fastq command to write all segments to FASTQ
        core                Core instruction. Imported by the rest.
        multiplex           Multiplex job file for each lane
        estimate            Prior estimate job file for each lane
        interleave          Interleaved job file for each lane
```
>

## IO API

`pheniqs-io-api.py` can manipulate a output format and splitting layout of Pheniqs configuration file.

>```shell
    usage: pheniqs-io-api.py [-h] -c PATH [-p PREFIX] [-f {fastq,sam,bam,cram}]
                             [--compression {uncompressed,gz,bgzf}]
                             [--compression-level {0,1,2,3,4,5,6,7,8,9}] [-l] [-s]
                             [-S] [-I PATH] [-O PATH] [--static] [--version]
                             [-v LEVEL]

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      -c PATH, --configuration PATH
                            Path to original configuration file.
      -p PREFIX, --prefix PREFIX
                            Prefix for generated output file names
      -f {fastq,sam,bam,cram}, --format {fastq,sam,bam,cram}
                            Output format
      --compression {uncompressed,gz,bgzf}
                            Output compression
      --compression-level {0,1,2,3,4,5,6,7,8,9}
                            Output compression level
      -l, --split-library   Library output routing
      -s, --split-segment   Segment output routing
      -S, --sense-input     sense input directive for pheniqs
      -I PATH, --base-input PATH
                            Base input URL directive for pheniqs
      -O PATH, --base-output PATH
                            Base output URL directive for pheniqs
      --static              Static configuratio output
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            logging verbosity level
```

## Prior API

`pheniqs-prior-api.py` can use an exiting configuration file and a Pheniqs demultiplex report for the same data to compile a new configuration file with adjusted priors.

>```shell
    usage: pheniqs-prior-api.py [-h] -c PATH [-r PATH] [-i [PATH [PATH ...]]]
                                [-I PATH] [-O PATH] [-s] [-p PREFIX] [--version]
                                [-v LEVEL]

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      -c PATH, --configuration PATH
                            Path to original configuration file.
      -r PATH, --report PATH
                            Path to report file.
      -i [PATH [PATH ...]], --input [PATH [PATH ...]]
                            Input directive for pheniqs
      -I PATH, --base-input PATH
                            Base input URL directive for pheniqs
      -O PATH, --base-output PATH
                            Base output URL directive for pheniqs
      -s, --sense-input     sense input directive for pheniqs
      -p PREFIX, --prefix PREFIX
                            Prefix for generated output file names
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            Logging verbosity level.
```
