#!/usr/bin/env Rscript

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti <lior.galanti@nyu.edu>
# Author: Alan Twaddle < alan.twaddle@nyu.edu >

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

# install.packages("ggplot2")
# install.packages("gridExtra")

library(grid)
library(ggplot2)
library(gridExtra)

source("theme.R")

diagram_width = 178
diagram_height = 60

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

benchmark_name_levels = c (
    "picard mdd split fastq",
    "picard mdd split cram",
    "fastq-multx mdd split fastq",
    "pheniqs mdd split fastq",
    "pheniqs pamld split fastq",
    "pheniqs mdd interleaved fastq",
    "pheniqs pamld interleaved fastq",
    "pheniqs mdd interleaved cram",
    "pheniqs pamld interleaved cram",
    "pheniqs mdd combined cram",
    "pheniqs pamld combined cram")

data = read.csv(data_filename, header=T, sep="\t")
data$name <- factor(data$description, 
    levels = benchmark_name_levels,
    labels = benchmark_name_labels
)
data$name <- with(data, reorder(name, quality))
quality_fill = scale_fill_manual (
    name = "quality",
    labels = c (
        "With",
        "Without"
    ),
    breaks = c (
        "0",
        "1"
    ),
    values = runtime_benchmark_color
)
flowcell_labeller = labeller (
    flowcell = c (
        "HCJFNBCXX" = "Dual Indexed 250bp Paired End",
        "HG7CVAFXX" = "Dual Indexed 75bp Paired End",
        "HK5NHBGXX" = "Dual Indexed 36bp Paired End"
    )
)
benchmark_plot = ggplot(data, aes(x = factor(name), y = seconds, fill=factor(quality), width=.75)) + 
benchmark_plot_theme + 
quality_fill +
ylab("Run time in minutes") +
scale_y_continuous(labels = function(x) round(x / 60, 0)) +
geom_bar(stat = "identity", position = "stack") +
facet_grid(. ~ flowcell, scales="free", labeller=flowcell_labeller) + coord_flip()

benchmark_grob <- ggplotGrob(benchmark_plot)
diagram <- arrangeGrob(benchmark_grob, ncol=1)
grid.draw(diagram)
if(grepl("\\.eps$", diagram_filename, perl = TRUE)) {
    ggsave(
        diagram_filename,
        diagram,
        scale = 1,
        width = diagram_width,
        height = diagram_height,
        units = diagram_units,
        device = cairo_ps,
        antialias = "subpixel")
}
if(grepl("\\.pdf$", diagram_filename, perl = TRUE)) {
    ggsave(
        diagram_filename,
        diagram,
        scale = 1,
        width = diagram_width,
        height = diagram_height,
        units = diagram_units,
        device = cairo_pdf)
}