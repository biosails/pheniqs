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

diagram_width = 86 * 3 * 0.6
diagram_height =  72 * 3  * 0.6

args = commandArgs(trailingOnly = TRUE)
diagram_filename = args[1]

variable_order = c (
    "FDR",
    "MR",
    "FP",
    "FN",
    "TP",
    "precision",
    "recall"
)

qc_order = c (
    "both",
    "fail",
    "pass"
)

tool_order = c (
    "deml",
    "mdd",
    "pamld"
)

variable_labeller = labeller (
    variable = c (
        "FDR" = "False discovery rate",
        "MR" = "Miss rate",
        "FP" = "False Positive",
        "FN" = "False Negative",
        "TP" = "True Positive",
        "precision" = "Precision",
        "recall" = "Recall"
    )
)

tool_color = scale_color_manual (
    name = "tool",
    labels = c (
        "deml",
        "mdd",
        "pamld"
    ),
    breaks = c (
        "deml",
        "mdd",
        "pamld"
    ),
    values = c (
        "deml" = alpha("#EB6841", 0.875),
        "mdd" = alpha("#151297", 0.875),
        "pamld" = alpha("#889725", 0.875)
    )
)
plot_measure <- function(data) {
    selected <- data
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$variable == 'FN'),]
    # selected <- selected[which(selected$qc == 'pass'),]
    # selected <- selected[which(selected$index == 0),]
    # selected <- selected[which(selected$value > 0),]
    # selected <- selected[which(selected$qc == 'fail'),]
    benchmark_plot <- ggplot(selected) +
    accuracy_plot_theme +
    facet_wrap(qc ~ variable, labeller = variable_labeller, scales="free") +
    # facet_grid(qc ~ variable, labeller = variable_labeller, scales="free") +
    # geom_line (
    #     data = selected,
    #     aes(
    #         x = rate,
    #         y = value,
    #         linetype = tool,
    #         colour = tool
    #     ),
    #     alpha = 0.5,
    #     size = 0.25
    # ) +
    geom_point (
        data = selected,
        aes(
            x = rate,
            y = value,
            colour = tool
        ),
        shape = 20,
        size = 0.5,
        alpha = 0.5,
        stroke = 0.75
    ) +
    # geom_jitter() +
    tool_color +
    # scale_linetype_manual (
    #     name = "tool",
    #     values = c (
    #         "pheniqs" = "solid",
    #         "deml" = "31"
    #     )
    # ) +
    # scale_y_log10() +
    # scale_x_log10() +
    guides(
        linetype = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position="top"
        )
    )

    return(benchmark_plot)
}

A5KVK_data = read.table('multiplex_R.csv', header=T, sep=",")
A5KVK_data$value = as.numeric(A5KVK_data$value)
# A5KVK_data$value = A5KVK_data$value + 0.000000000001
A5KVK_data$variable = factor(A5KVK_data$variable, levels = variable_order)
A5KVK_data$qc = factor(A5KVK_data$qc, levels = qc_order)
A5KVK_data$tool = factor(A5KVK_data$tool)
A5KVK_data$rate = as.numeric(A5KVK_data$rate)
A5KVK_plot <- plot_measure(A5KVK_data)
A5KVK_plot <- A5KVK_plot +
theme ( legend.position = "none" ) +
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
