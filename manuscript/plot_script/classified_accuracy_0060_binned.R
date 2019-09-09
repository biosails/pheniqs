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
maximum_error_rate = 0.0060
exclude_mdd = TRUE

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  bin = bin_name,
  ssid = experiment_id_name
)

rate_scale = scale_x_continuous (
    breaks = c (
        # 0.00010285707764030215,
        # 0.0005081994200929769,
        # 0.0010005303459722477,
        0.0010168814073261993, # not simulated
        # 0.0015038010220697069,
        0.0026127871615326784,
        0.005027553131321483,
        # 0.006922762955309007,
        0.009595885064725982,
        0.013481136183436082,
        0.01608628191772856,
        0.01927224171392631,
        0.023176701815888906,
        0.034072931974215626,
        0.04156162437196815,
        0.05084681116026185,
        0.0620952919053785,
        0.07516804510486658,
        0.10878556956710851,
        0.13098986579957753,
        0.15830773392739628
        # 0.19162814226990083
        # 0.23143485735751052
    ),
    labels = c (
        # 0.000103,
        # 0.0005082,
        # 0.001001,
        0.001017,
        # 0.001504,
        0.002613,
        0.005028,
        # 0.006923,
        0.009596,
        0.013481,
        0.016086,
        0.019272,
        0.023176,
        0.034073,
        0.041562,
        0.050847,
        0.062095,
        0.075168,
        0.108786,
        0.131990,
        0.158308
        # 0.191628
        # 0.231434
    )
)

plot_measure <- function(data) {
    selected <- data
    selected <- selected[which(selected$requested != 0),]
    selected <- selected[which(selected$rate < maximum_error_rate),]
    if(exclude_mdd) {
    selected <- selected[which(selected$tool != 'mdd'),]
    }
    # selected <- selected[which(selected$tool != 'deml'),]
    # selected <- selected[which(selected$tool != 'pamld_ap'),]
    # selected <- selected[which(selected$tool != 'pamld'),]
    # selected <- selected[which(selected$tool != 'pamld_u'),]
    # selected <- selected[which(selected$variable == 'FN'),]
    # selected <- selected[which(selected$variable == 'FP'),]
    # selected <- selected[which(selected$variable == 'TP'),]
    # selected <- selected[which(selected$variable == 'MR'),]
    # selected <- selected[which(selected$variable == 'FDR'),]
    # selected <- selected[which(selected$variable == 'fscore'),]

    selected <- selected[which(selected$variable != 'TP_FN'),]
    selected <- selected[which(selected$variable != 'TP_FP'),]
    selected <- selected[which(selected$variable != 'FN'),]
    selected <- selected[which(selected$variable != 'FP'),]
    selected <- selected[which(selected$variable != 'TP'),]
    selected <- selected[which(selected$variable != 'precision'),]
    selected <- selected[which(selected$variable != 'recall'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    # selected <- selected[which(selected$variable != 'fscore'),]

    benchmark_plot <- ggplot(selected) +
    pheniqs_plot_theme +
    theme(
      legend.position = "top",
      text = element_text(size = 10),
      axis.title.y = element_blank()
    ) +
    facet_wrap(
      bin ~ variable,
      labeller = accurecy_variable_labeller,
      scales = "free_y",
      strip.position = "left",
      ncol = 3
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
        size = 1,
        alpha = 0.325,
        stroke = 0.15
    ) +
    tool_color_scale +
    # scale_y_log10() +
    rate_scale
    return(benchmark_plot)
}

data = read.table(data_filename, header=T, sep=",")
data$value = as.numeric(data$value)
data$variable = factor(data$variable, levels = accurecy_variable_order)
data$tool = factor(data$tool)
data$bin = factor(data$bin)
data$rate = as.numeric(data$rate)
plot <- plot_measure(data)
plot <- plot +
# ggtitle( "Binned accuracy overview for classified reads only" ) +
xlab( "Substitution error rate" )

draw_diagram(plot, diagram_filename, diagram_width, diagram_height)
