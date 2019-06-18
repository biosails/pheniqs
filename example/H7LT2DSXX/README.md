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

# NovaSeq 6000
{:.page-title}

This is an example folder for processing a dual index paired end run performed by a NovaSeq 6000. The flowcell has 4 lanes each with 4 segments containing 151, 8, 8, and 151 cycles.

**All shell commands mentioned bellow are executed in this folder.**

The [181014_A00534_0024_AH7LT2DSXX](Illumia/181014_A00534_0024_AH7LT2DSXX) folder contains relevant metadata files from the run.

[core.json](core.json) is a core configuration file that summarized metadata extracted from the run and will be imported into the other configuration files. It was created by executing

```
illumina2pheniqs.py core illumina/181014_A00534_0024_AH7LT2DSXX > core.json
```

## Sample barcode decoding
The [sample](sample) folder contains configuration files and reports for estimating the sample barcode priors and decode them.

The [uniform](sample/uniform) folder contains configuration files for decoding the sample barcode with uniform priors.

The [prior](sample/prior) folder contains configuration files for decoding only the sample barcode with uniform priors and no output. Those are used for collecting statistics for the prior estimation and will generally execute much faster than full decoding. Reports for each run are also included and are used for estimating the priors.

The [adjusted](adjusted) folder contains configuration files for decoding the sample barcode with priors estimated from the repot. It is generated with the following commands

```
estimate_prior.py --report sample/prior/l01_sample_report.json --configuration sample/uniform/l01_sample.json > sample/adjusted/l01_sample.json
estimate_prior.py --report sample/prior/l02_sample_report.json --configuration sample/uniform/l02_sample.json > sample/adjusted/l02_sample.json
estimate_prior.py --report sample/prior/l03_sample_report.json --configuration sample/uniform/l03_sample.json > sample/adjusted/l03_sample.json
estimate_prior.py --report sample/prior/l04_sample_report.json --configuration sample/uniform/l04_sample.json > sample/adjusted/l04_sample.json
```
