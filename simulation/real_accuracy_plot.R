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
diagram_height =  72 * 1

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

variable_order = c (
    "FDR",
    "MR",
    "FP",
    "FN",
    "TP",
    "precision",
    "recall",
    "fscore"
)

qc_order = c (
    "both",
    "fail",
    "pass"
)

tool_order = c (
    "deml",
    "mdd",
    "pamld",
    "pamld_u",
    "pamld_ap"
)

variable_labeller = labeller (
    variable = c (
        "FDR" = "False discovery rate",
        "MR" = "Miss rate",
        "FP" = "False Positive",
        "FN" = "False Negative",
        "TP" = "True Positive",
        "precision" = "Precision",
        "recall" = "Recall",
        "fscore" = "F score"
    )
)

tool_color = scale_color_manual (
    name = "tool",
    labels = c (
        "deml",
        "mdd",
        "pamld",
        "pamld_u",
        "pamld_ap"
    ),
    breaks = c (
        "deml",
        "mdd",
        "pamld",
        "pamld_u",
        "pamld_ap"
    ),
    values = c (
        "deml" = alpha("#151297", 0.875),
        "mdd" = alpha("#9f0070", 0.875),
        "pamld" = alpha("#889725", 0.875),
        "pamld_u" = alpha("#f18407", 0.875),
        "pamld_ap" = alpha("#420d38", 0.875)
    )
)
plot_measure <- function(data) {
    selected <- data
    # selected <- selected[which(!(selected$qc == 'pass' & (selected$variable == 'FDR' | selected$variable == 'FP' | selected$variable == 'FN' | selected$variable == 'MR'))),]
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    selected <- selected[which(selected$rate < 0.06),]
    selected <- selected[which(selected$variable != 'FN'),]
    selected <- selected[which(selected$variable != 'FP'),]
    selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    selected <- selected[which(selected$qc == 'pass'),]
    # selected <- selected[which(selected$index == 0),]
    # selected <- selected[which(selected$value > 0),]
    # selected <- selected[which(selected$qc == 'fail'),]
    benchmark_plot <- ggplot(selected) +
    accuracy_plot_theme +
    facet_wrap(qc ~ variable, labeller = variable_labeller, scales="free") +
    # facet_grid(qc ~ variable, labeller = variable_labeller, scales="free") +
    geom_point (
        data = selected,
        aes(x = rate, y = value, colour = tool),
        shape = 21,
        size = 0.75,
        alpha = 0.75,
        stroke = 0.75
    ) +
    tool_color +
    # scale_x_log10() +
    # scale_y_log10() +
    # scale_x_log10 (
    #     breaks = c (
    #         0.00011424526022636,
    #         0.00052159516349246,
    #         0.00100105068181946,
    #         0.00104834072806070,
    #         0.00149408453025533,
    #         0.00260329363466613,
    #         0.00506003494781239,
    #         0.00704387579854628,
    #         0.00958068243750691,
    #         0.01148118043446911,
    #         0.01607627463062260,
    #         0.01613995231665040,
    #         0.01919226728146848,
    #         0.02323860965391178,
    #         0.05098662955880728,
    #         0.07541217255154616,
    #         0.10898342320002846,
    #         0.13086841380756978,
    #         0.15802835155326733,
    #         0.19174990026950645,
    #         0.19180936399101770,
    #         0.23159480764656637
    #     ),
    #     labels = c (
    #         "0.00011",
    #         "0.00052",
    #         "*0.00100",
    #         "0.00105",
    #         "0.00149",
    #         "0.00260",
    #         "0.00506",
    #         "0.00704",
    #         "0.00958",
    #         "0.01148",
    #         "0.01608",
    #         "0.01614",
    #         "0.01919",
    #         "0.02324",
    #         "0.05099",
    #         "0.07541",
    #         "0.10898",
    #         "0.13087",
    #         "0.15803",
    #         "0.19175",
    #         "0.19181",
    #         "0.23159"
    #     )
    # ) +
    scale_x_continuous (
        breaks = c (
            0.00010125944580606455,
            0.002184103418610344,
            0.0058997905077794125,
            0.009595387812854815,
            0.013468700997184169,
            0.016073157848319655,
            0.01925428239422401,
            0.023196602675082188,
            0.028042391655813122,
            0.03406639400032091,
            0.04152898980207327,
            0.05081449925444787,
            0.0620732865203883,
            0.07518373063317639
        ),
        labels = c (
            "0.00010",
            "0.00218",
            "0.00590",
            "0.00960",
            "0.01347",
            "0.01607",
            "0.01925",
            "0.02320",
            "0.02804",
            "0.03407",
            "0.04153",
            "0.05081",
            "0.06207",
            "0.07518"
        )
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
A5KVK_data$value = as.numeric(A5KVK_data$value)
A5KVK_data$variable = factor(A5KVK_data$variable, levels = variable_order)
A5KVK_data$qc = factor(A5KVK_data$qc, levels = qc_order)
A5KVK_data$tool = factor(A5KVK_data$tool)
A5KVK_data$rate = as.numeric(A5KVK_data$rate)
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
