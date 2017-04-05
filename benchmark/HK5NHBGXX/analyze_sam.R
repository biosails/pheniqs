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

library(grid)
library(ggplot2)
library(gridExtra)
library(extrafont)

source("../theme.R")

diagram_width = 178
diagram_height = 178

args = commandArgs(trailingOnly = TRUE)
sam_data_filename = args[1]
error_data_filename = args[2]
diagram_filename = args[3]

organism_names = c (
    "SC" = "S.cerevisiae",
    "NC" = "N.castellii",
    "CO" = "C.orthopsilosis",
    "SP" = "S.cerevisiae 2µ plasmid",
    "PX" = "PhiX",
    "BC" = "Bacteria",
    "PH" = "Phage",
    "VR" = "Virus",
    "UK" = "Unknown"
)
organism_labeller = labeller (
    OG = organism_names
)
organism_fill = scale_fill_manual (
    name = "Classification", 
    labels = organism_names, 
    breaks = c (
        "SC",
        "NC",
        "CO",
        "SP",
        "PX",
        "BC",
        "PH",
        "VR",
        "UK"
     ), 
     values = organism_color
)
organism_order = c ( 
    "UK",
    "VR",
    "PH",
    "BC",
    "PX",
    "SP",
    "CO",
    "NC",
    "SC"
)

# SAM plot

# CC DX WL WO OG AS MQ QS OBS OBQ BEE PBS PBD DBS DBD P CP REE RS RQ ID BS BD
# BS, DX, OG
# 22, 2,  5

# read and restructure error sam data
sam_data = read.csv(sam_data_filename, header=T, sep=" ")
clean_sam_data = sam_data[!(sam_data$BS %in% failed_libraries),]
classified = clean_sam_data[which(!grepl("=",clean_sam_data$BS)),]
classified$BS = factor(classified$BS)
classified$OG = factor(classified$OG, levels = organism_order)
classified$DX <- factor(classified$DX,
    levels = c (
        "0",
        "2",
        "1" 
    ),
    labels = c(
        "0" = "Consensus",
        "2" = expression(paste("P"[only])), 
        "1" = expression(paste("M"[only]))
    )
) 



# plot
difference_plot = ggplot(classified, aes(BS, fill=OG, width=0.75)) + 
geom_bar() +
sam_plot_theme +
organism_fill +
facet_grid(DX ~ ., labeller = label_parsed, scales="free_y") +
ggtitle("minimum distance vs phred-adjusted maximum likelihood decoding") +
xlab("Barcode") +
ylab("Read count") +
guides(fill = guide_legend(nrow=1)) +
scale_y_continuous()
difference_plot_grob <- ggplotGrob(difference_plot)

# ERROR plot
#
# read and restructure error data
error_data = read.csv(error_data_filename, header=T, sep="\t")
clean_error_data = error_data[!(error_data$base %in% failed_libraries | error_data$other %in% failed_libraries),]
clean_error_data$DX <- factor(clean_error_data$DX,
    levels = c (
        "0",
        "2",
        "1" 
    ),
    labels = c(
        "0" = "Consensus",
        "2" = expression(paste("P"[only], " vs Consensus")), 
        "1" = expression(paste("M"[only], " vs Consensus"))
    )
) 

self_error_data = clean_error_data[!(clean_error_data$self == 0),]
other_error_data = clean_error_data[!(clean_error_data$self == 1),]

error_plot = ggplot(other_error_data, aes(factor(base), R)) + 
error_plot_theme + 
geom_boxplot (
    lwd = 0.125,
    fatten = 1,
    width = 0.5,
    outlier.size = 0.375,
    color = box_color,
    fill = box_fill_color,
    outlier.color = outlier_color ) + 
geom_point(
    size = 0.625,
    data = self_error_data,
    color = self_point_color) +
facet_grid(DX ~ ., labeller = label_parsed) +
scale_y_log10() +
xlab("Barcode") +
ylab("Residual error")
error_plot_grob <- ggplotGrob(error_plot)

# align the two plots
maxWidth = grid::unit.pmax(error_plot_grob$widths, difference_plot_grob$widths)
error_plot_grob$widths <- as.list(maxWidth)
difference_plot$widths <- as.list(maxWidth)

diagram <- arrangeGrob(difference_plot_grob, error_plot_grob, ncol=1)
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