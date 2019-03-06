#!/usr/bin/env Rscript

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti < lior.galanti@nyu.edu >
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
# install.packages("extrafont")

library(grid)
library(ggplot2)
library(reshape2)
library(gridExtra)

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

source("theme.R")

diagram_width = 86 * 1
diagram_height =  72 * 1

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  rank = c (
    "real" = "",
    "noise" = "",
    "both" = ""
  ),
  qc = c (
      "pass" = "",
      "fail" = "",
      "both" = ""
  )
)

shift = 0.0000001

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$rate < maximum_rate),]
    selected <- selected[which(selected$rank == 'both'),]
    selected <- selected[which(selected$qc == 'pass'),]

    # selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    # selected <- selected[which(selected$tool != 'pamld_u'),]

    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]

    selected.melt <- melt(selected, id.vars=c('tool','rank','qc', 'TP', 'FP', 'FN', 'FDR', 'MR', 'precision', 'recall', 'fscore'))
    # selected.melt <- selected.melt[selected.melt$variable == "rate",]
    selected.melt$shifted_FDR = selected.melt$FDR + shift
    selected.melt$shifted_MR = selected.melt$MR + shift
    comparison_plot <- ggplot(selected.melt) +
    pheniqs_plot_theme +
    facet_wrap(~rank, labeller = accurecy_variable_labeller, scales="free") +
    geom_point (
      data = selected.melt,
      aes(
        x = shifted_MR,
        y = shifted_FDR,
        group = tool,
        color = tool,
        size = value
      ),
      shape = 20,
      alpha = 0.25,
      stroke = 0.5
    ) +
    xlab("Miss Rate + 1e-07") +
    ylab("False Discovery Rate + 1e-07") +
    scale_x_log10() +
    scale_y_log10() +
    tool_color_scale +
    scale_size_continuous (
      name = "Error Rate",
      breaks = c (
          0.00010285707764030215,
          # 0.0005081994200929769,
          # 0.0010005303459722477,
          # 0.0015038010220697069,
          # 0.0026127871615326784,
          0.005027553131321483,
          # 0.006922762955309007,
          0.009595885064725982,
          # 0.013481136183436082,
          0.016069265569637572,
          # 0.01608642715638242,
          # 0.01927224171392631,
          0.023176701815888906,
          0.05084681116026185,
          0.07516804510486658,
          0.10878556956710851,
          0.13098986579957753,
          0.15830773392739628,
          0.19162814226990083
      ),
      labels = c (
          0.0001029,
          # 0.0005082,
          # 0.001001,
          # 0.001504,
          # 0.002613,
          0.005028,
          # 0.006923,
          0.009596,
          # 0.01348,
          0.01607,
          # 0.01609,
          # 0.01927,
          0.02318,
          0.05085,
          0.07517,
          0.1088,
          0.131,
          0.1583,
          0.1916
      )
    )
    guides(
        size = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position="top"
        )
    )

    return(comparison_plot)
}

data = read.table(data_filename, header=T, sep=",")
plot <- plot_diagram(data)
plot <- plot +
ggtitle( "Comparing FDR and MR for Both Real and Noise reads" )

sheet <- ggplotGrob(plot)
diagram <- arrangeGrob(sheet, ncol=1)
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
        scale = 1.5,
        width = diagram_width,
        height = diagram_height,
        units = diagram_units,
        device = cairo_pdf)
}
