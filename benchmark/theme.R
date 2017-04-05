# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
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

library(extrafont)

diagram_units = "mm"
self_point_color = alpha("#333333", 0.75)
box_color = alpha("#3A3A3A", 0.75)
outlier_color = alpha("#3A3A3A", 0.25)
ticks_color = alpha("#3A3A3A", 1)
text_color = alpha("#3A3A3A", 1)
major_grid_color = alpha("#3A3A3A", 0.25)
axis_color = alpha("#3A3A3A", 0.75)
box_fill_color = alpha("#3A3A3A", 0.125)
font_family = "Monaco"
text_for_legend = element_text (
    size = rel(0.5),
    colour = text_color,
    family = font_family)
text_for_axis = element_text (
    size = rel(0.5),
    colour = text_color,
    family = font_family)
vertical_text_for_axis = element_text (
    angle = 90, 
    size = rel(0.5), 
    colour = text_color, 
    family = font_family)
text_for_title = element_text (
    size = rel(0.5), 
    colour = text_color, 
    family = font_family)

precision_recall_text_for_title = element_text (
    size = rel(0.375), 
    colour = text_color, 
    family = font_family)
precision_recall_text_for_legend = element_text (
    angle = 90, 
    size = rel(0.375),
    colour = text_color,
    family = font_family)
precision_recall_vertical_text_for_axis = element_text (
    angle = 90, 
    size = rel(0.375), 
    colour = text_color, 
    family = font_family)
precision_recall_text_for_axis = element_text (
    size = rel(0.375),
    colour = text_color,
    family = font_family,
    margin = margin(0,0,0,0))

runtime_benchmark_color = alpha ( c (
    "0" = "#384725",
    "1" = "#98A725"
), 0.875)
accuracy_decoder_color = c (
    "PAMLD" = alpha("#444444", 0.875),
    "MDD" = alpha("#333333", 0.875)
)
organism_color = c (
    "SC" = alpha( "#889725", 0.9),
    "NC" = alpha( "#98A725", 0.9),
    "CO" = alpha( "#A8B725", 0.9),
    "SP" = alpha( "#B8C725", 0.9),
    "PX" = alpha( "#FFAA07", 0.9),
    "BC" = alpha( "#6A4A3C", 0.9),
    "PH" = alpha( "#AA111D", 0.9),
    "VR" = alpha( "#EB6841", 0.9),
    "UK" = alpha( "#C94A65", 0.9)
)

tick_line = element_line (
    colour = ticks_color,
    size = rel(0.25),
    linetype = "solid" )
grid_line = element_line (
    colour = major_grid_color,
    size = rel(0.125),
    linetype = "solid" )
axis_line = element_line (
    colour = axis_color,
    size = rel(0.25),
    linetype = "solid" )
accuracy_plot_theme = theme (
    plot.title = precision_recall_text_for_title,
    legend.position = "left",
    legend.title = element_blank(),
    legend.text = precision_recall_text_for_legend,
    legend.text.align=-1.5,
    axis.line = axis_line,
    axis.ticks = tick_line,
    axis.title.x = precision_recall_text_for_title,
    axis.text.x = precision_recall_vertical_text_for_axis,
    axis.title.y = element_blank(),
    axis.text.y = precision_recall_text_for_axis,
    strip.text = precision_recall_text_for_title,
    strip.background = element_blank(),
    panel.grid.major = grid_line,
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.5, "char"),
    legend.key.height = unit(0.5, "char"),
    legend.key = element_rect(fill = NA),
    legend.margin = margin(0,0,0,0),
    legend.justification = c(0,0.5)
)
benchmark_plot_theme = theme (
    plot.title = element_blank(),
    legend.position="none",
    legend.title = element_blank(),
    legend.text = text_for_legend,
    axis.line = axis_line,
    axis.ticks = tick_line,
    axis.title.x = text_for_title,
    axis.text.x = vertical_text_for_axis,
    axis.title.y = element_blank(),
    axis.text.y = text_for_axis,
    strip.text = text_for_title,
    strip.background = element_blank(),
    panel.grid.major = grid_line,
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
)
sam_plot_theme = theme (
    plot.title = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = text_for_legend,
    axis.line = element_blank(),
    axis.ticks = tick_line,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = text_for_title,
    axis.text.y = text_for_axis,
    strip.text = text_for_title,
    strip.background = element_blank(),
    panel.grid.major = grid_line,
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.key.width = unit(0.375, "char"),
    legend.key.height = unit(0.5, "char"),
)
error_plot_theme = theme (
    plot.title = element_blank(),
    legend.position="none",
    legend.text = text_for_axis,
    axis.line = element_blank(),
    axis.ticks = tick_line,
    axis.title.x = text_for_title,
    axis.text.x = vertical_text_for_axis,
    axis.title.y = text_for_title,
    axis.text.y = text_for_axis,
    strip.text = text_for_title,
    strip.background = element_blank(),
    panel.grid.major = grid_line,
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
)
failed_libraries = c (
    # "TAAGGCGAAGGCTTAG",
    # "CGAGGCTGTACTCCTT",
    # "CAGAGAGGTACTCCTT",
    # "GTAGAGGATACTCCTT",
    # "TAAGGCGAAGAGGATA",

    "AAGAGGCAAGAGGATA",
    "AAGAGGCAAGGCTTAG",
    "AAGAGGCATATGCAGT",
    "AAGAGGCATCTTACGC",
    "CAGAGAGGAGGCTTAG",
    "CAGAGAGGTATGCAGT",
    "CAGAGAGGTCTTACGC",
    "CGAGGCTGAGAGGATA",
    "CGAGGCTGAGGCTTAG",
    "CGAGGCTGTATGCAGT",
    "CGTACTAGTCTACTCT",
    "CTCTCTACAGGCTTAG",
    "CTCTCTACTATGCAGT",
    "GCTACGCTAGAGGATA",
    "GCTACGCTAGGCTTAG",
    "GCTACGCTTATGCAGT",
    "GCTACGCTTCTTACGC",
    "GGACTCCTAGGCTTAG",
    "GGACTCCTTATGCAGT",
    "GTAGAGGAAGAGGATA",
    "GTAGAGGAAGGCTTAG",
    "GTAGAGGATATGCAGT",
    "GTAGAGGATCTTACGC",
    "TAAGGCGATACTCCTT",
    "TAAGGCGATCTACTCT",
    "TAGGCATGAGGCTTAG",
    "TAGGCATGTATGCAGT",
    "TCCTGAGCAGGCTTAG",
    "AGGCAGAATACTCCTT",
    "CGTACTAGTACTCCTT"
)
benchmark_name_labels = c (
    "picard fastq",
    "picard cram",
    "fastq-multx",
    "mdd split fastq",
    "pamld split fastq",
    "mdd interleaved fastq",
    "pamld interleaved fastq",
    "mdd interleaved cram",
    "pamld interleaved cram",
    "mdd combined cram",
    "pamld combined cram"
)
