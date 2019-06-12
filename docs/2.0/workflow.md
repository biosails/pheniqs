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

## Read anatomy design variants for Illumina platform

![paird end sequencing](/pheniqs/assets/img/paired_end_sequencing.png)

<a name="standard_illumina" />
Standard Illumina dual index, paird end sample demultiplexing
: This [vignette](illumina.html) will walk you through demultiplexing a standard Illumia high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
Classifying single index fluidigm with a single cellular tag
: This [vignette](fluidigm.html) will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).
