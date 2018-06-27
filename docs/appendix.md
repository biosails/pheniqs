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
        <li><a                  href="/pheniqs/">Home</a></li>
        <li><a                  href="/pheniqs/tutorial.html">Tutorial</a></li>
        <li><a                  href="/pheniqs/manual.html">Documentation</a></li>
        <li><a                  href="/pheniqs/glossary.html">Glossary</a></li>
        <li><a                  href="/pheniqs/workflow.html">Workflow</a></li>
        <li><a                  href="/pheniqs/cli.html">Command line interface</a></li>
        <li><a                  href="/pheniqs/building.html">Building</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>

#Appendix
{:.page-title}

##fastq-multx inefficient handling of dual indexing

[fastq-multx](https://github.com/brwnj/fastq-multx) handles dual indexing by applying a minimum distance decoder to the concatenated sequence of the two observed barcodes. While this may seem like a reasonable shortcut at first it is, in fact, inaccurate in some cases due to an unfortunate property of the the minimum distance decoding mismatch upper bound.

Consider an example of demultiplexing a **96** multiplex run with two **8bp** barcodes for a total of **16bp**.

When examining the minimum pairwise hamming distance, each of the two barcode sets can tolerate one error.
but when examining the pairwise hamming distance of a single barcode set composed of the the concatenated barcode, it can only tolerate one error.

Since fastq-multx examines the concatenated barcode it can only take one mismatch parameter and will not be able to handle the case where each of the two barcodes has one mismatch.
```
-d  Require a minimum distance of N between the best and next best
```

Setting the `-d` fastq-multx parameter to **1**, fastq-multx will declare a decoding error for the read.
Setting `-d` to **2** yields identical results because allowing **2** mismatches **anywhere** on the concatenated can yield ambiguous results.

To clarify, here is an example of the two observed barcode sequences where each had one mismatch and so the concatenated sequence had two:

```
Correct assignment      CTCTCTACAGAGGATA
observed sequence       CACTCTACAGAGGATT
observed quality        A/AAAAEAAA6AA/A/
mismatch                .*.............*
Error probability       1.50891e-07
Conditioned probability 0.00149889
```
**Notice the two mismatches are are on the lowest quality cycles of the base call**, which intuitively explains what went wrong. The probability of observing a random correctable barcode is **1.430511474e-06**
The conditioned probability is **1047** more likely than the strictest threshold. The probability of an incorrect assignment is **1.50891e-07**, which is **1** in **6,627,300**.

In a sample of **1,000,000** reads this occured **7622** times, or in **0.7622%** of the reads.
