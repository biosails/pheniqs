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
        <li><a                  href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a class="active"   href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>


# The Pheniqs Command Line
{:.page-title}

# Command line parameters
Command line parameters, if specified, override their corresponding values provided in the configuration file.

# zsh completion

If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you may wish to install the [bundled zsh command line completion script]({{ site.github.repository_url }}/blob/master/zsh/_pheniqs) for a more interactive command line experience.

# Global command line help

    pheniqs version 0.9.d55d65203560ae3719fcbb79ac657189e463ffed
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2017

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

    pheniqs version 0.9.d55d65203560ae3719fcbb79ac657189e463ffed
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2017

    Demultiplex and report quality control

    Usage : pheniqs demux [-h] [-V] [-D] -C PATH [-c FLOAT] [-f] [-q] [-n FLOAT]
                          [-m INT] [-i STRING] [-o STRING] [-l INT] [-d pamld|mdd|benchmark]
                          [-p CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-t INT]
                          [-T INT] [-L] [-B INT]

    Optional:
      -h, --help                  Show this help
      -V, --validate              Only validate configuration
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
      -L, --long                  Optimize threading for long read
      -B, --buffer INT            Records per resolution in feed buffer

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.

# Quality sub command

    pheniqs version 0.9.d55d65203560ae3719fcbb79ac657189e463ffed
    Lior Galanti < lior.galanti@nyu.edu >
    NYU Center for Genomics & Systems Biology 2017

    Report quality control

    Usage : pheniqs quality [-h] [-V] -i PATH [-f] [-L]
                            [-p CAPILLARY|LS454|ILLUMINA|SOLID|HELICOS|IONTORRENT|ONT|PACBIO] [-t INT]
                            [-T INT] [-B INT]

    Optional:
      -h, --help               Show this help
      -V, --validate           Only validate configuration
      -i, --input PATH         Path to input file
      -f, --filtered           Include filtered reads
      -L, --long               Optimize threading for long read
      -p, --platform STRING    Sequencing platform
      -t, --threads INT        IO thread pool size
      -T, --transforms INT     Number of transforming threads
      -B, --buffer INT         Records per resolution in feed buffer

    This program comes with ABSOLUTELY NO WARRANTY. This is free software,
    and you are welcome to redistribute it under certain conditions.