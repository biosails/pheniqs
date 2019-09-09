#!/usr/bin/env Rscript

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti < lior.galanti@nyu.edu >
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

library(reshape2)

accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  ssid = experiment_id_name
)

shift = 0.0000001

plot_diagram <- function(data) {
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

    # selected <- selected[which(selected$variable != 'FN'),]
    # selected <- selected[which(selected$variable != 'FP'),]
    # selected <- selected[which(selected$variable != 'TP'),]
    # selected <- selected[which(selected$variable != 'MR'),]
    # selected <- selected[which(selected$variable != 'FDR'),]
    selected.melt <- melt (
        selected,
        id.vars = c (
          'ssid',
          'expected',
          'requested',
          'tool',
          'classifier',
          'TP',
          'FP',
          'FN',
          'TP_FN',
          'TP_FP',
          'FDR',
          'MR',
          'precision',
          'recall',
          'fscore'
        )
    )
    selected.melt <- selected.melt[selected.melt$variable == "rate",]
    selected.melt$shifted_FDR = selected.melt$FDR + shift
    selected.melt$shifted_MR = selected.melt$MR + shift
    comparison_plot <- ggplot(selected.melt) +
    pheniqs_plot_theme +
    theme(
      legend.position = "left"
    ) +
    geom_point (
      data = selected.melt,
      aes(x = shifted_MR, y = shifted_FDR, group = tool, color = tool, size = value),
      shape = 20,
      alpha = 0.25,
      stroke = 0.5
    ) +
    xlab("Miss Rate + 1e-07") +
    ylab("False Discovery Rate + 1e-07") +
    scale_x_log10 (
      limits = c (shift, NA),
      breaks = c (
        0,
        0.00000001,
        0.0000001,
        0.000001,
        0.00001,
        0.0001,
        0.001,
        0.01,
        0.1,
        0.2
      ),
      labels = c (
        "0",
        "1e-08",
        "1e-07",
        "1e-06",
        "1e-05",
        "1e-04",
        "1e-03",
        "1e-02",
        "1e-01",
        "2e-01"
      )
    ) +
    scale_y_log10 (
      limits = c (0.000001, NA),
      breaks = c (
        0.000001,
        0.00001,
        0.00002,
        0.00003,
        0.00004,
        0.00006,
        0.00008,
        0.0001,
        0.0002,
        0.0003,
        0.0004,
        0.0006,
        0.0008,
        0.001,
        0.002,
        0.003,
        0.004,
        0.006,
        0.008,
        0.01,
        0.02
      )
    ) +
    tool_color_scale +
    scale_size_continuous (
      name = "Error Rate",
      breaks = c (
        0.00010285707764030215,
        # 0.0005081994200929769,
        # 0.0010005303459722477,
        # 0.0010168814073261993, # not simulated
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
        # 0.04156162437196815,
        0.05084681116026185,
        # 0.0620952919053785,
        0.07516804510486658
        # 0.10878556956710851,
        # 0.13098986579957753,
        # 0.15830773392739628
        # 0.19162814226990083
        # 0.23143485735751052
      ),
      labels = c (
        0.000103,
        # 0.0005082,
        # 0.001001,
        # 0.001017
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
        # 0.041562,
        0.050847,
        # 0.062095,
        0.075168
        # 0.108786,
        # 0.131990,
        # 0.158308
        # 0.191628
        # 0.231434
      )
    ) +
    guides(
        size = guide_legend (
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.position = "top"
        )
    )

    return(comparison_plot)
}

data = read.table(data_filename, header=T, sep=",")
plot <- plot_diagram(data)
plot <- plot +
ggtitle( "Classification false discovery and miss rates comparison" )

draw_diagram(plot, diagram_filename, diagram_width, diagram_height)
