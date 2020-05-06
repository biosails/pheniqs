---
layout: default
title: "Command Line Interface"
permalink: /cli/
id: cli
---

# The Pheniqs Command Line
{:.page-title}

# Command line parameters
Command line parameters, if specified, override their corresponding values provided in the configuration file.

# zsh completion
If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to install the bundled zsh command line completion script for a more interactive command line experience. It will interactively complete the command line arguments for you and makes learning the interface more intuitive. The zsh completion script, `_pheniqs`, is generated when building pheniqs with `make` but you can also generate it by invoking the corresponding `make` target `make _pheniqs` or with `shell.py` by executing `shell.py zsh configuration.json > _pheniqs`. The script should be placed in one of the folders referenced in your in your **fpath**.

# Global command line help

    pheniqs version 2.0.6
    For more information please visit http://biosails.github.io/pheniqs
    Lior Galanti < lior.galanti@nyu.edu > NYU Center for Genomics & Systems Biology 2018.

    Usage : pheniqs [-h] [--version] ACTION ...
    Optional :
      -h, --help    Show this help
      --version     Show program version


    Action :
      mux    Multiplex and Demultiplex annotated DNA sequence reads

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# `mux` sub command help

    pheniqs version 2.0.6-144
    Multiplex and Demultiplex annotated DNA sequence reads

    Usage : pheniqs mux [-h] [-i PATH]* [-o PATH]* [-c PATH] [-R PATH] [-I URL]
                        [-O URL] [-s] [-n] [-N] [-l INT] [-F fastq|sam|bam|cram] [-Z none|gz|bgzf]
                        [-L 0|1|2|3|4|5|6|7|8|9] [-T SEGMENT:START:END]*
                        [-P CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-q] [-V]
                        [-D] [-C] [-S] [-j] [-t INT] [--decoding-threads INT] [--htslib-threads INT]
                        [-B INT] [--precision INT]
    Optional :
      -h, --help                       Show this help
      -i, --input PATH                 Path to an input file. May be repeated.
      -o, --output PATH                Path to an output file. May be repeated.
      -c, --config PATH                Path to configuration file
      -R, --report PATH                Path to report file
      -I, --base-input URL             Base input url
      -O, --base-output URL            Base output url
      -s, --sense-input                Sense input segment layout
      -n, --no-output-npf              Filter outgoing QC failed reads
      -N, --no-input-npf               Filter incoming QC failed reads.
      -l, --leading INT                Leading read segment index
      -F, --format STRING              Defult output format
      -Z, --compression STRING         Defult output compression
      -L, --level STRING               Defult output compression level
      -T, --token SEGMENT:START:END    Output read token
      -P, --platform STRING            Sequencing platform
      -q, --quality                    Enable quality control
      -V, --validate                   Validate configuration file and emit a report
      -D, --distance                   Display pairwise barcode distance during validation
      -C, --compile                    Compiled JSON configuration file
      -S, --static                     Static configuration JSON file
      -j, --job                        Include a copy of the compiled job in the report
      -t, --threads INT                Thread pool size
      --decoding-threads INT           Number of parallel decoding threads
      --htslib-threads INT             Size of htslib thread pool size
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
      i.e. `pheniqs mux -i in_segment_1.fastq -i in_segment_2.fastq -o out_segment_1.fastq -o out_segment_2.fastq`

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Pheniqs tools
In the tool folder you will find several python scripts to assist with Pheniqs deployment and interfacing with existing tools.

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

* **basecall** sub command generates a shell script with a bcl2fastq command and an accompanying sample sheet for basecaling with illumina *bcl2fastq* without demupltiplexing.

* **core** generates a configuration that contains a single multiplex decoder directive in the global `decoder` directive for each lane in the sample sheet or just one if the integer `Lane` column is missing from the sample sheet.

* **multiplex** generates a multiplex job file for each lane with a uniform prior.

* **estimate** generates a prior estimation optimized job file for each lane.

* **interleave** generates an interleave job file for each lane. This can be used to pack the fastq files into a SAM format or just Interleaved fastq.

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

# JSON validation

JSON can be a little picky about syntax and a good JSON linter can make identifying offending syntax much easier. Plenty of tools for validating JSON syntax are out there but a simple good and readily available linter is available with the python programing language.

You will find a small [JSON linting python script]({{ site.github.repository_url }}/blob/master/tool/json_lint.py) in the tool directory that is somewhat customized for the pheniqs config file.

For **python 2** use:
>```shell
    python -c "import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')"
```

or for **python 3**:
>```shell
    python3 -c "import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))"
```

You may alternatively set it up as an alias in your shell's profile by adding to your `.zshrc` or `.bashrc`:
>```shell
    alias jsl="python -c \"import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')\""
```

or

>```shell
    alias jsl="python3 -c \"import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))\""
```

and than invoke it simply by feeding it a JSON file on standard input:
>```shell
    cat configuration.json|jsl
```

This will print an easy to read tabulated JSON to standard output and assist you with resolving any syntactical JSON violations.
