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


#Appendix

##zsh completion

If you use [zsh](https://en.wikipedia.org/wiki/Z_shell) you will find a [command line completion](zsh/_pheniqs) script bundled with the code. It will interactively complete the command line arguments for you and makes learning the interface more intuitive. It should be placed or symlinked in a folder that is in your **fpath**.

##The JSON configuration file
JSON can be a little picky about syntax and a good JSON linter can make identifying offending syntax much easier. Plenty of tools for validating JSON syntax are out there but a simple good and readily available linter is available with the python programing language.

For **python 2** use:

```python -c "import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')"```

or for **python 3**:

```python3 -c "import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))"```

You may alternatively set it up as an alias in your shell's profile by adding to your `.zshrc` or `.bashrc`:

```alias jsl="python -c \"import json,sys; print json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4).encode('utf8')\""```

or 

```alias jsl="python3 -c \"import json,sys; print(json.dumps(json.load(sys.stdin),sort_keys=True,ensure_ascii=False,indent=4))\""```

and than invoke it simply by feeding it a JSON file on standard input:

```cat configuration.json|jsl```

This will print an easy to read tabulated JSON to standard output and assist you with resolving any syntactical JSON violations.

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
