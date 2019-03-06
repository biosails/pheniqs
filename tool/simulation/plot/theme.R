# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti <lior.galanti@nyu.edu>

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

library(extrafont)
library(Cairo)

maximum_rate = 0.06
diagram_units = "mm"
self_point_color = alpha("#333333", 0.75)
box_color = alpha("#3A3A3A", 0.75)
outlier_color = alpha("#3A3A3A", 0.25)
ticks_color = alpha("#3A3A3A", 1)
text_color = alpha("#3A3A3A", 1)
major_grid_color = alpha("#3A3A3A", 0.25)
axis_color = alpha("#3A3A3A", 0.75)
box_fill_color = alpha("#3A3A3A", 0.125)
font_family = "Monaco"

plot_title_text = element_text (
  size = rel(0.5),
  hjust = 0.5,
  colour = text_color,
  family = font_family )

title_text = element_text (
  size = rel(0.375),
  hjust = 0.5,
  colour = text_color,
  family = font_family )

axis_text = element_text (
  size = rel(0.375),
  colour = text_color,
  family = font_family,
  margin = margin(0,0,0,0) )

legend_text = element_text (
  angle = 0,
  size = rel(0.375),
  colour = text_color,
  family = font_family )

vertical_axis_text = element_text (
  angle = 90,
  size = rel(0.375),
  colour = text_color,
  family = font_family)

tick_line = element_line (
  colour = ticks_color,
  size = rel(0.25),
  linetype = "solid" )

grid_line = element_line (
  colour = major_grid_color,
  size = rel(0.125),
  linetype = "solid" )

axis_line = element_line (
  colour = axis_color,
  size = rel(0.25),
  linetype = "solid" )

pheniqs_plot_theme = theme (
  plot.title = plot_title_text,
  legend.position = "right",
  legend.title = legend_text,
  legend.text = legend_text,
  # legend.text.align=-1.5,
  axis.line = axis_line,
  axis.ticks = tick_line,
  axis.title.x = title_text,
  axis.text.x = vertical_axis_text,
  axis.title.y = vertical_axis_text,
  axis.text.y = axis_text,
  strip.text = title_text,
  strip.background = element_blank(),
  panel.grid.major = grid_line,
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key.width = unit(0.5, "char"),
  legend.key.height = unit(0.5, "char"),
  legend.key = element_rect(fill = NA),
  legend.margin = margin(0,0,0,0),
  legend.justification = c(0,0.5)
)

tool_order = c (
  "pamld",
  "pamld_u",
  "pamld_ap",
  "deml",
  "mdd"
)
tool_name = c (
  "pamld" = "PAMLD Estimated",
  "pamld_u" = "PAMLD Uniform",
  "pamld_ap" = "PAMLD Accurate",
  "deml" = "deML",
  "mdd" = "MDD"
)
tool_color = c (
  "pamld" = alpha("#889725", 0.875),
  "pamld_u" = alpha("#f18407", 0.875),
  "pamld_ap" = alpha("#420d38", 0.875),
  "deml" = alpha("#151297", 0.875),
  "mdd" = alpha("#9f0070", 0.875)
)
tool_linetype = c (
  "pamld" = "solid",
  "pamld_u" = "solid",
  "pamld_ap" = "solid",
  "deml" = "solid",
  "mdd" = "solid"
)
tool_color_scale = scale_color_manual (
  name = "Tool",
  breaks = tool_order,
  labels = tool_name,
  values = tool_color
)
tool_linetype_scale = scale_linetype_manual (
    name = "tool",
    values = tool_linetype
)
accurecy_variable_order = c (
  "FDR",
  "MR",
  "FP",
  "FN",
  "TP",
  "precision",
  "recall",
  "fscore"
)
accurecy_variable_name = c (
  "FDR" = "False Discovery Rate",
  "MR" = "Miss Rate",
  "FP" = "False Positive",
  "FN" = "False Negative",
  "TP" = "True Positive",
  "precision" = "Precision",
  "recall" = "Recall",
  "fscore" = "F score"
)
accurecy_rank_name = c (
  "real" = "Real reads only",
  "noise" = "Noise reads only",
  "both" = "Both Real and Noise reads"
)
quality_control_order = c (
  "both",
  "fail",
  "pass"
)
quality_control_name = c (
  "pass" = "Passed Quality Control",
  "fail" = "Failed Quality Control",
  "both" = "Both Passed and Failed Quality Control"
)
accurecy_variable_labeller = labeller (
  tool = tool_name,
  variable = accurecy_variable_name,
  rank = accurecy_rank_name,
  qc = quality_control_name
)

# rate_scale = scale_x_continuous (
rate_scale = scale_x_log10 (
    breaks = c (
        0.00010285707764030215,
        # 0.0005081994200929769,
        0.0010005303459722477,
        # 0.0015038010220697069,
        # 0.0026127871615326784,
        0.005027553131321483,
        # 0.006922762955309007,
        0.009595885064725982,
        # 0.013481136183436082,
        0.016069265569637572,
        # 0.01608642715638242,
        # 0.01927224171392631,
        0.023176701815888906,
        0.05084681116026185,
        0.07516804510486658,
        0.10878556956710851,
        0.13098986579957753,
        0.15830773392739628
        # 0.19162814226990083
    ),
    labels = c (
        0.0001029,
        # 0.0005082,
        0.001001,
        # 0.001504,
        # 0.002613,
        0.005028,
        # 0.006923,
        0.009596,
        # 0.01348,
        0.01607,
        # 0.01609,
        # 0.01927,
        0.02318,
        0.05085,
        0.07517,
        0.1088,
        0.131,
        0.1583
        # 0.1916
    )
)
