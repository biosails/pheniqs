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

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$rate < maximum_error_rate),]
    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "none",
      text = element_text(size = 10)
    ) +
    geom_col(
      data = selected,
      aes(x = ssid, y = error, fill = index),
      alpha = 0.75,
      size = 0.25
    )
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$error = as.numeric(data$error)
data$ssid = factor(data$ssid, labels = experiment_id_name, levels = experiment_id_order)
plot <- plot_diagram(data)
plot <- plot +
ggtitle( "Prior estimation error" ) +
xlab( "Substitution error rate" ) +
ylab( "Error" )

draw_diagram(plot, diagram_filename, diagram_width, diagram_height)
