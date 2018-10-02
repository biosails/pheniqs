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
        <li><a                  href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/best_practices.html">Best Practice</a></li>
        <li><a                  href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a class="active"   href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# The Pheniqs Command Line
{:.page-title}

# Command line parameters
Command line parameters, if specified, override their corresponding values provided in the configuration file.

# zsh completion
If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to install the bundled zsh command line completion script for a more interactive command line experience. It will interactively complete the command line arguments for you and makes learning the interface more intuitive. The zsh completion script, `_pheniqs`, is generated when building pheniqs with `make` but you can also generate it by invoking the corresponding `make` target `make _pheniqs` or with `pheniqs-tools.py` by executing `pheniqs-tools.py zsh configuration.json > _pheniqs`. The script should be placed in one of the folders referenced in your in your **fpath**.

# Global command line help

    pheniqs version 2.0.4
    For more information please visit http://biosails.github.io/pheniqs
    Lior Galanti < lior.galanti@nyu.edu > NYU Center for Genomics & Systems Biology 2018.

    Usage : pheniqs [-h] [--version] ACTION ...
    Optional :
      -h, --help    Show this help
      --version     Show program version


    Action :
      demux    Multiplex and Demultiplex annotated DNA sequence reads

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Demux sub command help

    pheniqs version 2.0.4
    Multiplex and Demultiplex annotated DNA sequence reads

    Usage : pheniqs demux [-h] [-i PATH]* [-o PATH]* [-c PATH] [-R PATH] [-I URL]
                          [-O URL] [-s] [-f] [-l INT] [-F fastq|sam|bam|cram] [-Z none|gz]
                          [-T SEGMENT:START:END]*
                          [-P CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-q] [-V]
                          [-D] [-C] [-S] [-j] [-t INT] [-B INT] [--precision INT]
    Optional :
      -h, --help                       Show this help
      -i, --input PATH                 Path to an input file. May be repeated.
      -o, --output PATH                Path to an output file. May be repeated.
      -c, --config PATH                Path to configuration file
      -R, --report PATH                Path to report file
      -I, --base-input URL             Base input url
      -O, --base-output URL            Base output url
      -s, --sense-input                Sense input segment layout
      -f, --filtered                   Include reads not passing vendor QC in output
      -l, --leading INT                Leading read segment index
      -F, --format STRING              Defult output format
      -Z, --compression STRING         Defult output compression
      -T, --token SEGMENT:START:END    Output read token
      -P, --platform STRING            Sequencing platform
      -q, --quality                    Enable quality control
      -V, --validate                   Validate configuration file and emit a report
      -D, --distance                   Display pairwise barcode distance during validation
      -C, --compile                    Compiled JSON configuration file
      -S, --static                     Static configuration JSON file
      -j, --job                        Include a copy of the compiled job in the report
      -t, --threads INT                Thread pool size
      -B, --buffer INT                 Feed buffer capacity
      --precision INT                  Output floating point precision

      -i/--input defaults to /dev/stdin with inputing layout sensing.
      -o/--output default to /dev/stdout with SAM format.
      -I/--base-input and -O/--base-output default to the working directory.
      -V/--validate, -C/--compile and -S/--static disable job excution and only emit information.
      -s/--sense-input will guess input layout by examining the first few reads of each input file.
      -S/--static emits a static configuration file with all imports resolved.
      -C/--compile emits a compiled configuration file ready for execution with implicit attributes resolved.
      -i/--input and -o/--output can be repeated to provide multiple paths,
      i.e. `pheniqs demux -i in_segment_1.fastq -i in_segment_2.fastq -o out_segment_1.fastq -o out_segment_2.fastq`

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Pheniqs tools
In the tool folder you will find several python scripts to assist with Pheniqs deployment and interfacing with existing tools.

## `illumina2pheniqs.py`

Generate pheniqs configuration files or a bcl2fastq command from an illumina run directory. This tool parses that `RunInfo.xml`, `RunParameters.xml` and `SampleSheet.csv` files in the directory. The `Data` section of the `SampleSheet.csv` must either have all records declare a `Lane` or none.

    usage: illumina2pheniqs.py [-h] [--version] [-v LEVEL] ACTION ...

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
        bcl2fastq           bcl2fastq command to write all segments to FASTQ
        core                Single decoder directive for each lane
        interleave          Interleave both template and index segments to SAM
        demultiplex         Demultiplex a single lane

The `core` sub command generates a configuration that contains a single multiplex decoder directive in the global `decoder` directive for each lane
in the sample sheet or just one if the integer `Lane` column is missing from the sample sheet.
Those declarations are not used in the `multiplex` directive but are merely made available for other configurations to import, where the exact output layout and possible additional decoders can be declared.

The `demultiplex` sub command will declare an inline multiplex decoder for a single lane, specified with the `-l/--lane-number` command line parameter.
not specifying a lane will create a multiplex decoder with all records for sample sheets that do not declare a `Lane`, like the MiSeq.

The `interleave` sub command generates a configuration file that will interleave all FASTQ files containing the segments of the read into a single stream. By default this will emit the segments to stdout in SAM format but you may explicitly specify an output file path in any other format.

The `bcl2fastq` sub command generates a shell command for executing `bcl2fastq 2.x` that will disable all post processing and demultiplexing done with bcl2fastq and only convert the data in the bcl files into a FASTQ file for each segment.

# JSON validation

JSON can be a little picky about syntax and a good JSON linter can make identifying offending syntax much easier. Plenty of tools for validating JSON syntax are out there but a simple good and readily available linter is available with the python programing language.

You will find a small [JSON linting python script]({{ site.github.repository_url }}/blob/master/tool/json_lint.py) in the tool directory that is somewhat customized for the pheniqs config file.

For **python 2** use:

    python -c "import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')"

or for **python 3**:

    python3 -c "import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))"

You may alternatively set it up as an alias in your shell's profile by adding to your `.zshrc` or `.bashrc`:

    alias jsl="python -c \"import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')\""

or

    alias jsl="python3 -c \"import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))\""

and than invoke it simply by feeding it a JSON file on standard input:

    cat configuration.json|jsl

This will print an easy to read tabulated JSON to standard output and assist you with resolving any syntactical JSON violations.
