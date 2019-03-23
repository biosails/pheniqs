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

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

source("theme.R")

diagram_width = 86 * 1
diagram_height =  72 * 1

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  bin = bin_name,
  ssid = experiment_id_name
)

plot_diagram <- function(data) {
    selected <- data

    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "top",
      text = element_text(size = 8)
    ) +
    geom_bar(
      data = selected,
      stat = "identity",
      position = position_dodge(),
      aes(x = order, y = count, fill = bin),
      # fill = alpha("#5C5151", 0.875),
      alpha = 0.5,
      size = 0.25
    )
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$bin = factor(data$bin, labels = bin_name)
plot <- plot_diagram(data)
plot <- plot +
ggtitle( "Library density distribution" ) +
xlab( "Library" ) +
ylab( "Read Count" )

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
