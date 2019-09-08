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

diagram_width = 170
diagram_height =  170
maximum_error_rate = 0.0550

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$rate < maximum_error_rate),]
    selected <- selected[which(selected$density > 0),]
    selected <- selected[which(selected$density < 20000000),]
    benchmark_plot <- ggplot(selected) +
    facet_wrap(~ssid, ncol = 2) +
    pheniqs_plot_theme +
    theme(
        legend.position = "none",
        text = element_text(size = 12)
    ) +
    geom_col(
      data = selected,
      position = position_dodge(),
      aes(x = quality, y = density, fill = quality),
      alpha = 0.7,
      size = 0.25
    ) +
    phred_quality_fill_scale
    # scale_y_log10() +
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$quality = factor(data$quality)
data$density = as.numeric(data$density)
data$ssid = factor(data$ssid, labels = experiment_id_name, levels = experiment_id_order)
plot <- plot_diagram(data)
plot <- plot +
# ggtitle( "Quality distribution for barcode cycles" ) +
xlab( "Quality" ) +
ylab( "Read ount" )

draw_diagram(plot, diagram_filename, diagram_width, diagram_height)
