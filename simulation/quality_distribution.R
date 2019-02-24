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
# install.packages("extrafont")

library(grid)
library(ggplot2)
library(gridExtra)

source("theme.R")

diagram_width = 86 * 5
diagram_height =  72 * 5

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

plot_measure <- function(data) {
    selected <- data
    benchmark_plot <- ggplot(selected) +
    facet_wrap(vars(rate)) +
    geom_bar(
      data = selected,
      stat = "identity",
      position = position_dodge(),
      aes(
          x = quality,
          y = density
      ),
      alpha = 0.5,
      size = 0.25,
      show.legend = FALSE
    ) +
    guides(
        linetype = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position="top"
        )
    )
    return(benchmark_plot)
}

A5KVK_data = read.table(data_filename, header=T, sep=",")
A5KVK_data$quality = as.numeric(A5KVK_data$quality)
A5KVK_data$density = as.numeric(A5KVK_data$density)
A5KVK_data$rate = factor(A5KVK_data$rate)
A5KVK_plot <- plot_measure(A5KVK_data)
A5KVK_plot <- A5KVK_plot +
# theme ( legend.position = "none" ) +
ggtitle( "100 libraries / 2x7bp barcode" ) +
xlab( "Expected error rate" )
A5KVK_grob <- ggplotGrob(A5KVK_plot)

diagram <- arrangeGrob(A5KVK_grob, ncol=1)
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
