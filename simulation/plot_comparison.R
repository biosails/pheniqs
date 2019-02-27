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

source("theme.R")

diagram_width = 86 * 2
diagram_height =  72 * 1

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$rank != 'noise'),]
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$tool != 'deml'),]
    selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    selected <- selected[which(selected$tool != 'pamld_u'),]
    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    selected <- selected[which(selected$qc == 'pass'),]
    # selected <- selected[selected$rate > 0.002,]

    selected.melt <- melt(selected, id.vars=c('tool','rank','qc', 'TP', 'FP', 'FN', 'FDR', 'MR'))
    selected.melt <- selected.melt[selected.melt$variable == "rate",]

    comparison_plot <- ggplot(selected.melt) +
    facet_wrap(~rank, labeller = accurecy_variable_labeller, scales="free") +
    comparison_plot_theme +
    # geom_line (
    #   data = selected.melt,
    #   aes(
    #     x = FP,
    #     y = TP,
    #     group = tool,
    #     color = tool,
    #     size = value
    #   ),
    #   alpha = 0.5,
    #   size = 0.25,
    #   show.legend = FALSE
    # ) +
    geom_point (
      data = selected.melt,
      aes(
        x = FN,
        y = TP,
        group = tool,
        color = tool,
        size = value
      ),
      shape = 20,
      alpha = 0.3,
      stroke = 0.5
    ) +
    # scale_x_continuous(limits = c(0, NA)) +
    # scale_x_log10() +
    tool_color +
    guides(
        linetype = guide_legend (
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
# theme ( legend.position = "none" ) +
# xlab( "Expected Nucleotide Error Rate" )
ggtitle( "100 libraries / 2x7bp barcode" )
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
        scale = 1,
        width = diagram_width,
        height = diagram_height,
        units = diagram_units,
        device = cairo_pdf)
}
