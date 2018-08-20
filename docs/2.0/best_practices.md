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

# Tips
Pheniqs is a generic, fast, and accurate sequence demultiplexer, apt to 
demultiplex almost every conceivable barcoding strategy in use today.
That said, there are some common conventions we recommend for the best possible 
performance, regardless of the barcoding strategy employed. Below are a few tips 
that we recommend that we believe are best for your analysis.

## Reduce all read segments into CRAM BEFORE Analysis
[CRAM](glossary.html#htslib) is an efficient format for storing sequencing information. 
Since Pheniqs is I/O bound, reading from only one file is significantly faster.

- Automatic input detection
- 
