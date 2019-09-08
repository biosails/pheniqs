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

source("core.R")
library("lubridate", warn.conflicts = FALSE)

diagram_width = 85
diagram_height =  170

duration_scale = scale_y_continuous (
  name = "Duration",
  breaks = c(
      0,
      3600 * 2,
      3600 * 4,
      3600 * 8,
      3600 * 16,
      3600 * 32,
      3600 * 64,
      287460
  ),
  labels = c(
      "00:00:00",
      "02:00:00",
      "04:00:00",
      "08:00:00",
      "16:00:00",
      "32:00:00",
      "64:00:00",
      "79:51:49"
  )
)
configuration_order = c (
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7"
    # "8",
    # "9"
)
configuration_name = c (
    "deml fastq SI to SS fastq",
    "Pheniqs fastq SI to null",
    "Pheniqs fastq SI to II bam",
    "Pheniqs fastq SI to SI fastq",
    "Pheniqs fastq SI to SS fastq",
    "Pheniqs fastq SI to II cram",
    "Pheniqs mdd fastq SI to II bam"
)
configuration_color = c (
  "1" = alpha("#C00606", 1),
  "2" = alpha("#3852BB", 1),
  "3" = alpha("#3852BB", 1),
  "4" = alpha("#3852BB", 1),
  "5" = alpha("#3852BB", 1),
  "6" = alpha("#3852BB", 1),
  "7" = alpha("#3852BB", 1)
)
configuration_scale = scale_fill_manual (
  name = "Configuration",
  values = configuration_color,
  breaks = configuration_order
)

plot_speed_diagram <- function(data) {
    selected <- data
    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "none",
      axis.text.x = diagonal_axis_text,
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.title.y = element_blank(),
      text = element_text(size = 10)
    ) +
    geom_col(
      data = selected,
      aes(x = index, y = time, fill = tool),
      alpha = 0.75,
      size = 0.25
    ) +
    # coord_flip() +
    duration_scale +
    configuration_scale
    return(benchmark_plot)
}

plot_memory_diagram <- function(data) {
    selected <- data
    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.title.y = element_blank(),
      text = element_text(size = 10)
    ) +
    geom_col(
      data = selected,
      aes(x = index, y = memory, fill = tool),
      alpha = 0.75,
      size = 0.25
    ) +
    # coord_flip() +
    configuration_scale
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data <- data[which(data$tool != 'bcl2fastq'),]
# data <- data[which(data$tool != 'deml'),]

data$tool = factor(data$index)
data$index = factor(data$index, labels = configuration_name, levels = configuration_order)
data$time = lubridate::seconds(lubridate::hms(data$duration))

speed_plot <- plot_speed_diagram(data)
speed_plot <- speed_plot +
ggtitle( "Speed" ) +
ylab( "Time" ) +
xlab( "Configuration" )
speed_plot_grob <- ggplotGrob(speed_plot)

memory_plot <- plot_memory_diagram(data)
memory_plot <- memory_plot +
ggtitle( "Memory" ) +
ylab( "Megabytes" ) +
xlab( "Configuration" )
memory_plot_grob <- ggplotGrob(memory_plot)

# align the two plots
maxWidth = grid::unit.pmax(speed_plot_grob$widths, memory_plot_grob$widths)
speed_plot_grob$widths <- as.list(maxWidth)
memory_plot_grob$widths <- as.list(maxWidth)
diagram <- arrangeGrob(memory_plot_grob, speed_plot_grob, ncol=1)

# maxHeight = grid::unit.pmax(speed_plot_grob$heights, memory_plot_grob$heights)
# speed_plot_grob$heights <- as.list(maxHeight)
# memory_plot_grob$heights <- as.list(maxHeight)
# diagram <- arrangeGrob(speed_plot_grob, memory_plot_grob, ncol=2)

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
