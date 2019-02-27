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

diagram_width = 86 * 3
diagram_height =  72 * 2

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

plot_measure <- function(data) {
    selected <- data
    selected <- selected[which(selected$rate < 0.1),]
    selected <- selected[which(selected$rank == 'both'),]
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$rate > 0.00104 | selected$rate < 0.0014),]
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    # selected <- selected[which(selected$tool != 'pamld_u'),]
    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    selected <- selected[which(selected$qc == 'pass'),]
    benchmark_plot <- ggplot(selected) +
    accuracy_plot_theme +
    facet_wrap(qc ~ variable, labeller = accurecy_variable_labeller, scales="free") +
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
    geom_point (
        data = selected,
        aes(x = rate, y = value, colour = tool),
        shape = 21,
        size = 0.75,
        alpha = 0.5,
        stroke = 0.5
    ) +
    tool_linetype +
    tool_color +
    scale_x_log10() +
    # scale_y_log10() +
    # scale_x_continuous (
    #     breaks = c (
    #         0.00010125944580606455,
    #         0.002184103418610344,
    #         0.0058997905077794125,
    #         0.009595387812854815,
    #         0.013468700997184169,
    #         0.016073157848319655,
    #         0.01925428239422401,
    #         0.023196602675082188,
    #         0.028042391655813122,
    #         0.03406639400032091,
    #         0.04152898980207327,
    #         0.05081449925444787,
    #         0.0620732865203883,
    #         0.07518373063317639
    #     ),
    #     labels = c (
    #         "0.00010",
    #         "0.00218",
    #         "0.00590",
    #         "0.00960",
    #         "0.01347",
    #         "0.01607",
    #         "0.01925",
    #         "0.02320",
    #         "0.02804",
    #         "0.03407",
    #         "0.04153",
    #         "0.05081",
    #         "0.06207",
    #         "0.07518"
    #     )
    # ) +
    guides(
        linetype = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position="top"
        )
    )
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$value = as.numeric(data$value)
data$variable = factor(data$variable, levels = accurecy_variable_order)
data$qc = factor(data$qc, levels = accurecy_qc_order)
data$tool = factor(data$tool)
data$rate = as.numeric(data$rate)
plot <- plot_measure(data)
plot <- plot +
# theme ( legend.position = "none" ) +
ggtitle( "100 libraries / 2x7bp barcode" ) +
xlab( "Expected error rate" )
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
