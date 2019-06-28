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
