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

    usage: pheniqs-illumina-api.py [-h] [--version] [-v LEVEL] ACTION ...

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            logging verbosity level

    pipeline operations:
      Generate pheniqs configuration files or a bcl2fastq command from an
      illumina run directory. This tool parses RunInfo.xml, RunParameters.xml
      and SampleSheet.csv files.

      ACTION
        basecall            bcl2fastq command to write all segments to FASTQ
        core                Core instruction. Imported by the rest.
        multiplex           Multiplex job file for each lane
        estimate            Prior estimate job file for each lane
        interleave          Interleaved job file for each lane

## IO API

`pheniqs-io-api.py` can manipulate a output format and splitting layout of Pheniqs configuration file.

    usage: pheniqs-io-api.py [-h] -c PATH [-F {fastq,sam,bam,cram}]
                             [--compression {uncompressed,gz,bgzf}]
                             [--compression-level {0,1,2,3,4,5,6,7,8,9}] [-L] [-S]
                             [--base-input PATH] [--base-output PATH] [-s]
                             [--static] [-p PREFIX] [--version] [-v LEVEL]

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      -c PATH, --configuration PATH
                            Path to original pheniqs configuration file.
      -F {fastq,sam,bam,cram}, --format {fastq,sam,bam,cram}
                            Output format
      --compression {uncompressed,gz,bgzf}
                            Output compression
      --compression-level {0,1,2,3,4,5,6,7,8,9}
                            Output compression level
      -L, --split-library   Library output routing
      -S, --split-segment   Segment output routing
      --base-input PATH     Forwarded to pheniqs -I/--base-input parameter.
      --base-output PATH    Forwarded to pheniqs -O/--base-output parameter.
      -s, --sense-input     sense input directive for pheniqs
      --static              Static configuration output
      -p PREFIX, --prefix PREFIX
                            Prefix for generated output file names
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            logging verbosity level

## Prior API

`pheniqs-prior-api.py` can use an exiting configuration file and a Pheniqs demultiplex report for the same data to compile a new configuration file with adjusted priors.

    usage: pheniqs-prior-api.py [-h] -c PATH [-r PATH] [-p PREFIX] [-i PATH] [-s]
                                [--base-input PATH] [--base-output PATH]
                                [--version] [-v LEVEL]

    Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology

    optional arguments:
      -h, --help            show this help message and exit
      -c PATH, --configuration PATH
                            Path to original pheniqs configuration file.
      -r PATH, --report PATH
                            Path to pheniqs prior estimation report file.
      -p PREFIX, --prefix PREFIX
                            Generated file names prefix.
      -i PATH, --input PATH
                            Forwarded to pheniqs -i/--input parameter.
      -s, --sense-input     Forwarded to pheniqs -s/--sense-input parameter.
      --base-input PATH     Forwarded to pheniqs -I/--base-input parameter.
      --base-output PATH    Forwarded to pheniqs -O/--base-output parameter.
      --version             show program's version number and exit
      -v LEVEL, --verbosity LEVEL
                            logging verbosity level
