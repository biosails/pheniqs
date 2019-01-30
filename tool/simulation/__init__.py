#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2018  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti <lior.galanti@nyu.edu>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Structure of a SAM record
# -----------------------------------------------------------------------------
# 0  QNAME   string     Query template NAME
# 1  FLAG    int        bitwise FLAG
#    0x1     template having multiple segments in sequencing
#    0x2     each segment properly aligned according to the aligner
#    0x4     segment unmapped
#    0x8     next segment in the template unmapped
#    0x10    SEQ being reverse complemented
#    0x20    SEQ of the next segment in the template being reverse complemented
#    0x40    the first segment in the template
#    0x80    the last segment in the template
#    0x100   secondary alignment
#    0x200   not passing filters, such as platform/vendor quality controls
#    0x400   PCR or optical duplicate
#    0x800   supplementary alignment
# 2  RNAME   string     Reference sequence NAME
# 3  POS     int        1-based leftmost mapping POSition
# 4  MAPQ    int        MAPping Quality
# 5  CIGAR   string     CIGAR string
# 6  RNEXT   string     Reference name of the mate/next read
# 7  PNEXT   int        Position of the mate/next read
# 8  TLEN    int        observed Template LENgth
# 9  SEQ     string     segment SEQuence
# 10 QUAL    string     Phred QUALity+33

from simulation.barcode import SimulateBarcode
from simulation.substitution import SimulateSubstitution
from simulation.deml import ToDeML
from simulation.demultiplex import PamldDemultiplex
from simulation.demultiplex import PamldAccuratePriorDemultiplex
from simulation.demultiplex import PamldUniformDemultiplex
from simulation.demultiplex import MddDemultiplex
from simulation.demultiplex import DemlDemultiplex
from simulation.analyze import Analyze
from simulation.prior import SensePrior, AdjustPrior
from simulation.summarize import Summarize
# from simulation.collect import Collect
