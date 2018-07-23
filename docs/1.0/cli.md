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
        <li><a                  href="/pheniqs/1.0/">Home</a></li>
        <li><a                  href="/pheniqs/1.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/1.0/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/1.0/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/1.0/workflow.html">Workflow</a></li>
        <li><a class="active"   href="/pheniqs/1.0/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/1.0/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>


# The Pheniqs Command Line
{:.page-title}

# Command line parameters
Command line parameters, if specified, override their corresponding values provided in the configuration file.

# zsh completion

If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to install the [bundled zsh command line completion script]({{ site.github.repository_url }}/blob/master/zsh/_pheniqs) for a more interactive command line experience. It will interactively complete the command line arguments for you and makes learning the interface more intuitive. It should be placed or symlinked in a folder that is in your **fpath**.

# Global command line help

    pheniqs version 0.9.d55d65203560ae3719fcbb79ac657189e463ffed
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2018

    Usage : pheniqs [-h] [--version] ACTION ...

    Optional:
      -h, --help    Show this help
      --version     Show program version

    available action
      demux      Demultiplex and report quality control
      quality    Report quality control

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Demux sub command help

    pheniqs 1.0.2068ec583e9ad5ab75877cd47ee977565c945d60 I’ll build my own theme park
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2018

    Demultiplex and report quality control

    Usage : pheniqs demux [-h] [-V] [-j] [-D] -C PATH [-c FLOAT] [-f] [-q]
                          [-n FLOAT] [-m INT] [-i STRING] [-o STRING] [-l INT] [-d pamld|mdd|benchmark]
                          [-p CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-t INT]
                          [-T INT] [-L] [-B INT]

    Optional:
      -h, --help                  Show this help
      -V, --validate              Only validate configuration
      -j, --lint                  Only lint configuration file
      -D, --distance              Display pairwise barcode distance
      -C, --config PATH           Path to configuration file
      -c, --confidence FLOAT      Decoding confidence threshold
      -f, --filtered              Include filtered reads
      -q, --quality               Disable quality control
      -n, --noise FLOAT           Noise prior probability
      -m, --mask INT              Phred masking threshold
      -i, --base-input STRING     Base input path
      -o, --base-output STRING    Base output path
      -l, --leading INT           Leading read segment
      -d, --decoder STRING        Barcode decoder
      -p, --platform STRING       Sequencing platform
      -t, --threads INT           IO thread pool size
      -T, --transforms INT        Number of transforming threads
      -B, --buffer INT            Records per resolution in feed buffer

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Quality sub command

    pheniqs 1.0.2068ec583e9ad5ab75877cd47ee977565c945d60 I’ll build my own theme park
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2018

    Report quality control

    Usage : pheniqs quality [-h] [-V] -i PATH [-f] [-L]
                            [-p CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-t INT]
                            [-T INT] [-B INT]

    Optional:
      -h, --help               Show this help
      -V, --validate           Only validate configuration
      -i, --input PATH         Path to input file
      -f, --filtered           Include filtered reads
      -p, --platform STRING    Sequencing platform
      -t, --threads INT        IO thread pool size
      -T, --transforms INT     Number of transforming threads
      -B, --buffer INT         Records per resolution in feed buffer

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# JSON validation

JSON can be a little picky about syntax and a good JSON linter can make identifying offending syntax much easier. Plenty of tools for validating JSON syntax are out there but a simple good and readily available linter is available with the python programing language.

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
