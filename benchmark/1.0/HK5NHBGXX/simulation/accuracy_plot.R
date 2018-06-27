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

source("../../theme.R")

diagram_width = 86
diagram_height =  72

args = commandArgs(trailingOnly = TRUE)
diagram_filename = args[1]

variable_order = c (
    "precision",
    "recall",
    "f_score"
)

variable_labeller = labeller (
    variable = c (
        "precision" = "False discovery rate",
        "recall" = "Miss rate",
        "f_score" = "F-score"
    )
)
plot_measure <- function(data, nature='natural', confidence_min=NULL, confidence_max=NULL) {
    selected <- data
    selected <- selected[which(selected$nature == nature),]
    if(!is.null(confidence_min)) {
        selected <- selected[which(selected$confidence >= confidence_min),]
    }
    if(!is.null(confidence_max)) {
        selected <- selected[which(selected$confidence <= confidence_max),]
    }
    benchmark_plot <- ggplot(selected) + 
    accuracy_plot_theme + 
    facet_grid(variable ~ ., labeller = variable_labeller, scales="free_y") +
    geom_line (
        data = selected[selected$decoder == "PAMLD",],
        aes (
            x = confidence,
            y = value,
            linetype = decoder
        ),
        alpha = 0.5,
        size = 0.25
    ) +
    geom_point (
        data = selected[selected$decoder == "PAMLD",],
        aes(
            x = confidence,
            y = value
        ),
        shape = 21,
        size = 1,
        alpha = 0.5,
        stroke = 0.5
    ) + 
    scale_linetype_manual (
        name = "decoder",
        values = c (
            "PAMLD" = "solid", 
            "MDD" = "31"
        )
    ) +
    scale_colour_manual (
        values = accuracy_decoder_color
    ) +
    scale_x_log10(
        breaks = c (
            0.000001,
            0.000010,
            0.000100,
            0.001000,
            0.010000,
            0.020000,
            0.050000,
            0.100000,
            0.200000,
            0.300000,
            1.000000
        )
        ,labels = c (
            "0.0001",
            "0.001",
            "0.01",
            "0.1",
            "1",
            "2",
            "5",
            "10",
            "20",
            "30",
            "100"
        )
    ) +
    guides(
        linetype = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position="top"
        )
    ) +
    geom_hline (
        data = selected[selected$decoder == "MDD",],
        aes(yintercept = value, color = decoder),
        size = 0.25,
        alpha = 0.75,
        linetype="31"
    )


    return(benchmark_plot)
}

SK5NHBGXY_data = read.table('SK5NHBGXY_measure.csv', header=T, sep="\t")
SK5NHBGXY_data$variable = factor(SK5NHBGXY_data$variable, levels = variable_order)
SK5NHBGXY_data$confidence = as.numeric(SK5NHBGXY_data$confidence)
SK5NHBGXY_plot <- plot_measure(SK5NHBGXY_data, 'natural', 0.0001, 0.1)
SK5NHBGXY_plot <- SK5NHBGXY_plot + 
theme ( legend.position = "none" ) + 
ggtitle( "10 libraries / 6bp barcode" ) + 
xlab( "% Permissible error" )
SK5NHBGXY_grob <- ggplotGrob(SK5NHBGXY_plot)

SK5NHBGXX_data = read.table('SK5NHBGXX_measure.csv', header=T, sep="\t")
SK5NHBGXX_data$variable = factor(SK5NHBGXX_data$variable, levels = variable_order)
SK5NHBGXX_data$confidence = as.numeric(SK5NHBGXX_data$confidence)
SK5NHBGXX_plot <- plot_measure(SK5NHBGXX_data, 'natural', 0.00001, 0.1)
SK5NHBGXX_plot <- SK5NHBGXX_plot + 
theme ( strip.text = element_blank(), legend.position = "none") + 
ggtitle( "96 libraries / dual 8bp barcode" ) + 
xlab( "% Permissible error" )
SK5NHBGXX_grob <- ggplotGrob(SK5NHBGXX_plot)

# align the two plots
maxHeight = grid::unit.pmax(SK5NHBGXY_grob$heights, SK5NHBGXX_grob$heights)
SK5NHBGXY_grob$heights <- as.list(maxHeight)
SK5NHBGXX_grob$heights <- as.list(maxHeight)

diagram <- arrangeGrob(SK5NHBGXX_grob, SK5NHBGXY_grob, ncol=2)
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