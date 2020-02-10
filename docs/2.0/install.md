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
        <li><a class="active"   href="/pheniqs/2.0/">Home</a></li>
        <li><a                  href="/pheniqs/2.0/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/2.0/build.html">Vignettes</a></li>
        <li><a class="active"   href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a                  href="/pheniqs/2.0/cli.html">CLI</a></li>
        <li><a                  href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

# Installation
{:.page-title}

Pheniqs is distributed as precompiled binaries, which may be installed with package managers as described below. To build it from scratch, please look <a href="/pheniqs/2.0/build.html">here</a>.

## Installing with homebrew on MacOS
**coming soon**

## Installing with conda

### One time setup - Install Miniconda
The easiest way to do this is to head on over to [anaconda](https://conda.io/miniconda.html) and select the correct distribution for python3.

### One time setup - Configure your channels
Many groups contribute software to conda. Each of these groups corresponds to a different channel. Bioconda is a well known channel for distributing bioinformatics software. It depends on conda-forge, another group for distributing more general software, including R and python packages.

>```shell
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

### Install Pheniqs
Simply install pheniqs using the conda package manager

>```shell
# Installs pheniqs in your 'global' conda install
conda install pheniqs
# Installs pheniqs in an isolated environment - recommended
conda create -n pheniqs pheniqs
```

If you want to live on the bleeding edge, you can also install pheniqs from our anaconda channel, nyuad-cgsb.

>```shell
conda install -c nyuad-cgsb pheniqs/latest
```
