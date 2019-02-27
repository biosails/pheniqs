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
text_for_legend = element_text (
    size = rel(0.5),
    colour = text_color,
    family = font_family)
text_for_axis = element_text (
    size = rel(0.5),
    colour = text_color,
    family = font_family)
vertical_text_for_axis = element_text (
    angle = 90,
    size = rel(0.5),
    colour = text_color,
    family = font_family)
text_for_title = element_text (
    size = rel(0.5),
    colour = text_color,
    family = font_family)

precision_recall_text_for_title = element_text (
    size = rel(0.375),
    colour = text_color,
    family = font_family)
precision_recall_text_for_legend = element_text (
    angle = 0,
    size = rel(0.375),
    colour = text_color,
    family = font_family)
precision_recall_vertical_text_for_axis = element_text (
    angle = 90,
    size = rel(0.375),
    colour = text_color,
    family = font_family)
precision_recall_text_for_axis = element_text (
    size = rel(0.375),
    colour = text_color,
    family = font_family,
    margin = margin(0,0,0,0))
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
accuracy_plot_theme = theme (
    plot.title = precision_recall_text_for_title,
    legend.position = "left",
    legend.title = element_blank(),
    legend.text = precision_recall_text_for_legend,
    # legend.text.align=-1.5,
    axis.line = axis_line,
    axis.ticks = tick_line,
    axis.title.x = precision_recall_text_for_title,
    axis.text.x = precision_recall_vertical_text_for_axis,
    axis.title.y = element_blank(),
    axis.text.y = precision_recall_text_for_axis,
    strip.text = precision_recall_text_for_title,
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
comparison_plot_theme = theme (
    plot.title = precision_recall_text_for_title,
    legend.position = "left",
    legend.title = element_blank(),
    legend.text = precision_recall_text_for_legend,
    # legend.text.align=-1.5,
    axis.line = axis_line,
    axis.ticks = tick_line,
    axis.title.x = precision_recall_text_for_title,
    axis.text.x = precision_recall_vertical_text_for_axis,
    axis.title.y = precision_recall_text_for_title,
    # axis.title.y = element_blank(),
    axis.text.y = precision_recall_text_for_axis,
    strip.text = precision_recall_text_for_title,
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
accurecy_variable_labeller = labeller (
    variable = accurecy_variable_name,
    rank = accurecy_rank_name,
    qc = c (
        "pass" = "Passed Quality Control",
        "fail" = "Failed Quality Control",
        "both" = "Both Passed and Failed Quality Control"
    )
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
accurecy_qc_order = c (
    "both",
    "fail",
    "pass"
)
accurecy_tool_order = c (
    "deml",
    "mdd",
    "pamld",
    "pamld_u",
    "pamld_ap"
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
tool_linetype = scale_linetype_manual (
    name = "tool",
    values = c (
      "deml" = "solid",
      "mdd" = "solid",
      "pamld" = "solid",
      "pamld_u" = "31",
      "pamld_ap" = "31"
    )
)
