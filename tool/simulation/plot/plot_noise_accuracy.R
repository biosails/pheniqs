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

diagram_width = 86 * 2
diagram_height =  72 * 3

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  rank = accurecy_rank_name,
  qc = c (
      "pass" = "",
      "fail" = "",
      "both" = ""
  )
)

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$rate < maximum_erro_rate),]
    selected <- selected[which(selected$qc == 'fail'),]
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    # selected <- selected[which(selected$tool != 'pamld_u'),]
    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      axis.title.y = element_blank()
    ) +
    facet_wrap(
      qc ~ variable,
      labeller = accurecy_variable_labeller,
      scales="free",
      ncol=2
    ) +
    geom_line (
        data = selected,
        aes(
            x = rate,
            y = value,
            linetype = tool,
            colour = tool
        ),
        alpha = 0.5,
        size = 0.25,
        show.legend = FALSE
    ) +
    tool_linetype_scale +
    geom_point (
        data = selected,
        aes(
            x = rate,
            y = value,
            colour = tool
        ),
        shape = 21,
        size = 1.25,
        alpha = 0.325,
        stroke = 0.25
    ) +
    tool_color_scale +
    rate_scale +
    guides(
        linetype = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position = "top"
        )
    )

    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$value = as.numeric(data$value)
data$variable = factor(data$variable, levels = accurecy_variable_order)
data$qc = factor(data$qc, levels = quality_control_order)
data$tool = factor(data$tool)
data$rate = as.numeric(data$rate)
plot <- plot_diagram(data)
plot <- plot +
# theme ( legend.position = "none" ) +
ggtitle( "Accuracy Statistics for Noise Reads Only" ) +
xlab( "Expected Nucleotide Error Rate" )
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
