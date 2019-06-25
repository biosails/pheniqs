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
<!-- <p>&#8636;&#8637;</p> -->
<table class="diagram">
  <tr>
    <td class="description" >Single index</td>
    <td>
      <div class="read" id="single_index">
        <div class="binding_primer p5">P5</div>
        <div class="sequencing_primer">SP1</div>
        <div class="insert">insert</div>
        <div class="sequencing_primer">SP2</div>
        <div class="sample_barcode i7">i7</div>
        <div class="binding_primer p7">P7</div>
        <div class="clear"></div>
      </div>
    </td>
  </tr>
  <tr>
    <td class="description" >Paired-end dual index</td>
    <td>
      <div class="read" id="single_index">
        <div class="binding_primer p5">P5</div>
        <div class="sample_barcode i5">i5</div>
        <div class="sequencing_primer">SP1</div>
        <div class="insert">insert</div>
        <div class="sequencing_primer">SP2</div>
        <div class="sample_barcode i7">i7</div>
        <div class="binding_primer p7">P7</div>
        <div class="clear"></div>
      </div>
    </td>
  </tr>
  <tr>
    <td class="description" >Paired-end dual index with UMI</td>
    <td>
      <div class="read" id="single_index">
        <div class="binding_primer p5">P5</div>
        <div class="sample_barcode i5">i5</div>
        <div class="sequencing_primer">SP1</div>
        <div class="insert">insert</div>
        <div class="sequencing_primer">SP2</div>
        <div class="sample_barcode i7">i7</div>
        <div class="molecular_barcode">UMI</div>
        <div class="binding_primer p7">P7</div>
        <div class="clear"></div>
      </div>
    </td>
  </tr>
  <tr>
    <td class="description" >Paired-end dual index with inline UMI</td>
    <td>
      <div class="read" id="single_index">
        <div class="binding_primer p5">P5</div>
        <div class="sample_barcode i5">i5</div>
        <div class="sequencing_primer">SP1</div>
        <div class="molecular_barcode">UMI</div>
        <div class="insert">insert</div>
        <div class="sequencing_primer">SP2</div>
        <div class="sample_barcode i7">i7</div>
        <div class="binding_primer p7">P7</div>
        <div class="clear"></div>
      </div>
    </td>
  </tr>
  <tr>
    <td class="description" >Paired-end dual index, inline UMI and barcode</td>
    <td>
      <div class="read" id="single_index">
        <div class="binding_primer p5">P5</div>
        <div class="sample_barcode i5">i5</div>
        <div class="sequencing_primer">SP1</div>
        <div class="molecular_barcode">UMI</div>
        <div class="custom_barcode">BC</div>
        <div class="insert">insert</div>
        <div class="sequencing_primer">SP2</div>
        <div class="sample_barcode i7">i7</div>
        <div class="binding_primer p7">P7</div>
        <div class="clear"></div>
      </div>
    </td>
  </tr>
</table>

<table class="legend">
  <tr>
    <td class="label" ><span class="binding_primer p5">P5/P7</span></td>
    <td class="definition" >Platform specific flow cell binding sequences</td>
  </tr>
  <tr>
    <td class="label" ><span class="sequencing_primer">SP1/SP2</span></td>
    <td class="definition" >Sequencing primer binding sites (common for all libraries)</td>
  </tr>
  <tr>
    <td class="label" ><span class="sample_barcode">i7/i5</span></td>
    <td class="definition" >Library specific sample indexes</td>
  </tr>
  <tr>
    <td class="label" ><span class="molecular_barcode">UMI</span></td>
    <td class="definition" >Unique molecular index (barcode tag for individual molecules)</td>
  </tr>
  <tr>
    <td class="label" ><span class="custom_barcode">BC</span></td>
    <td class="definition" >User-defined barcode (unique per sample, single cell, etc.)</td>
  </tr>
  <tr>
    <td class="label" ><span class="insert">insert</span></td>
    <td class="definition" >Target DNA or cDNA fragment (library-specific)</td>
  </tr>
</table>

<a name="novaseq_6000" />
Illumina NovaSeq6000 with `pheniqs-illumina-api`
: The [NovaSeq 6000 vignette](novaseq.html) will walk you through demultiplexing a dual indexed paired end NovaSeq 6000 run with the `pheniqs-illumina-api` python interface. It will show you how to generate configuration files from the Illumina run folder, estimate the sample barcode priors and demultiplex the run. It will apply to almost every standard sample multiplex Illumina run with minor changes.

<a name="standard_illumina" />
Standard Illumina dual index, paird end sample demultiplexing
: This [vignette](illumina.html) will walk you through demultiplexing a standard Illumia high throughput sequencing run with paired end dual index samples multiplexed using the standard Illumina i5 and i7 index protocol with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).

<a name="fluidigm" />
Classifying single index fluidigm with a single cellular tag
: This [vignette](fluidigm.html) will walk you through a single index fluidigm sequencing run with the [PAMLD decoder](glossary.html#phred_adjusted_maximum_likelihood_decoding).
