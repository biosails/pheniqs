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

source("theme.R")

args = commandArgs(trailingOnly = TRUE)
data_filename = args[1]
diagram_filename = args[2]

diagram_width = 86 * 2
diagram_height =  72 * 3

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  ssid = experiment_id_name
)

plot_diagram <- function(data) {
    selected <- data
    selected <- selected[which(selected$requested != 0),]
    selected <- selected[which(selected$rate < maximum_error_rate),]
    selected <- selected[which(selected$tool != 'mdd'),]
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    # selected <- selected[which(selected$tool != 'pamld_u'),]
    # selected <- selected[which(selected$variable != 'TP_FP'),]
    selected <- selected[which(selected$variable != 'TP_FN'),]
    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    selected <- selected[which(selected$variable != 'recall'),]
    selected <- selected[which(selected$variable != 'precision'),]
    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "top",
      axis.title.y = element_blank()
    ) +
    facet_wrap(
      ~ variable,
      labeller = accurecy_variable_labeller,
      strip.position = "left",
      scales="free",
      ncol=2
    ) +
    geom_line (
        data = selected,
        aes(x = rate, y = value, colour = tool),
        alpha = 0.5,
        size = 0.25
    ) +
    geom_point (
        data = selected,
        aes(x = rate, y = value, colour = tool),
        shape = 21,
        size = 1.25,
        alpha = 0.325,
        stroke = 0.25
    ) +
    tool_color_scale +
    rate_scale

    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$value = as.numeric(data$value)
data$variable = factor(data$variable, levels = accurecy_variable_order)
data$tool = factor(data$tool)
data$rate = as.numeric(data$rate)
plot <- plot_diagram(data)
plot <- plot +
ggtitle( "Classification accuracy for unclassified reads only" ) +
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
