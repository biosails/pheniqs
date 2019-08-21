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

# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("extrafont")
# install.packages("lubridate")

library(extrafont)
library(Cairo)
library(grid)
library(ggplot2)
library(gridExtra)

data_filename = NULL
diagram_filename = NULL
maximum_error_rate = 0.5
diagram_width = 86
diagram_height =  72

args = commandArgs(trailingOnly = TRUE)
if(length(args) > 0) { data_filename = args[1] }
if(length(args) > 1) { diagram_filename = args[2] }
if(length(args) > 2) { maximum_error_rate = as.numeric(args[3]) }

diagram_units = "mm"
font_family = "Monaco"
text_color = alpha("#3A2323", 1)
axis_color = alpha("#3A2323", 0.75)
ticks_color = alpha("#3A2323", 1)
major_grid_color = alpha("#3A2323", 0.25)

plot_title_text = element_text (
  size = rel(0.5),
  hjust = 0.5,
  colour = text_color,
  family = font_family
)
title_text = element_text (
  size = rel(0.375),
  hjust = 0.5,
  colour = text_color,
  family = font_family
)
axis_text = element_text (
  size = rel(0.375),
  colour = text_color,
  family = font_family,
  margin = margin(0,0,0,0)
)
legend_text = element_text (
  size = rel(0.375),
  colour = text_color,
  family = font_family,
  margin = margin(0,0,0,0)
)
legend_title_text = element_text (
  size = rel(0.375),
  colour = text_color,
  family = font_family,
  margin = margin(0,0,0,0, "char")
)
vertical_axis_text = element_text (
  angle = 90,
  size = rel(0.375),
  colour = text_color,
  family = font_family
)
diagonal_axis_text = element_text (
  angle = 30,
  hjust = 1,
  size = rel(0.375),
  colour = text_color,
  family = font_family
)
tick_line = element_line (
  colour = ticks_color,
  size = rel(0.25),
  linetype = "solid"
)
grid_line = element_line (
  colour = major_grid_color,
  size = rel(0.125),
  linetype = "solid"
)
axis_line = element_line (
  colour = axis_color,
  size = rel(0.25),
  linetype = "solid"
)
tool_order = c (
  "pamld",
  "pamld_u",
  "pamld_ap",
  "deml",
  "mdd"
)
tool_name = c (
  "pamld" = "PAMLD estimated",
  "pamld_u" = "PAMLD uniform",
  "pamld_ap" = "PAMLD true",
  "deml" = "deML",
  "mdd" = "MDD"
)
tool_color = c (
  "pamld" = alpha("#5F9A10", 1),
  "pamld_u" = alpha("#A080C7", 1),
  "pamld_ap" = alpha("#3852BB", 1),
  "deml" = alpha("#C00606", 1),
  "mdd" = alpha("#5C5151", 1)
)
tool_color_scale = scale_color_manual (
  name = "Tool",
  breaks = tool_order,
  labels = tool_name,
  values = tool_color,
)
accurecy_variable_order = c (
  "FDR",
  "MR",
  "fscore",
  "precision",
  "recall",
  "TP",
  "FP",
  "FN",
  "TP_FP",
  "TP_FN"
)
accurecy_variable_name = c (
  "FDR" = "False Discovery Rate",
  "MR" = "Miss Rate",
  "FP" = "False Positive",
  "FN" = "False Negative",
  "TP" = "True Positive",
  "precision" = "Precision",
  "recall" = "Recall",
  "fscore" = "F score",
  "TP_FP" = "True Positive + False Positive",
  "TP_FN" = "True Positive + False Negative"
)
bin_order = c (
  "0",
  "1",
  "2",
  "3",
  "4"
)
bin_color = c (
  "0" = alpha("#5F9A10", 1),
  "1" = alpha("#A080C7", 1),
  "2" = alpha("#3852BB", 1),
  "3" = alpha("#C00606", 1),
  "4" = alpha("#5C5151", 1)
)
bin_name = c (
  "0" = "% < 0.001",
  "1" = "0.001 < % < 0.003",
  "2" = "0.003 < % < 0.01",
  "3" = "0.01 < % < 0.03",
  "4" = "0.03 < %"
)
# bin_name = c (
#   "0" = "% < 0.001 / 23",
#   "1" = "0.001 < % < 0.003 / 44",
#   "2" = "0.003 < % < 0.01 / 5",
#   "3" = "0.01 < % < 0.03 / 27",
#   "4" = "0.03 < % / 1"
# )
bin_fill_scale = scale_fill_manual (
  name = "Bin",
  breaks = bin_order,
  labels = bin_name,
  values = bin_color,
  aesthetics = "fill"
)
experiment_id_order = c (
  "37e41038-310c-44e1-a1ed-8979b0f73709",
  "082e6a59-27fe-47ed-ae2c-89da8ceb31ee",
  "ce333fb9-5e8d-4b14-a5ee-7101afe2cfa6",
  "b3b1b9af-9a6c-43dd-b5ab-be35bbee9c7f",
  "40a08990-33b3-41f5-b061-2561bc84ca96",
  "8c689128-8483-4ce5-9b67-8cb160550497",
  "b7e78bf9-30a3-4955-9469-0bbd77310476",
  "9442e6ea-8033-44b2-b0f3-91470348a7f0",
  "156d68c8-abcd-439a-a387-71b3708652ff",
  "7d4caa7a-5b9c-4b47-af6a-0217ce344f6a",
  "a362969d-fe30-497a-a2ec-42bc0119f398",
  "62bdb0f5-affe-4b31-8b27-f5e96899f37f",
  "b37d56c6-488a-4121-abf6-8b8a2bf82668",
  "424c4ff7-5cc6-4032-bd54-70b492bbcbb5",
  "a5b2c471-f1bc-49b7-9650-19bcd78e73ae",
  "58020739-79fd-4cab-8f4e-e611a6cfba15",
  "ffad3de5-49d4-4c8d-8be9-6db1be42dbda",
  "92ffff7e-8c09-447a-96b1-9bcdc2a0da1c",
  "a11ae384-a142-4f33-a129-eda212675e2a",
  "44eb6f63-6f19-47d2-acce-67f2e59395bf",
  "1a0d8c5d-3a1b-48ec-894a-12888c735955",
  "3c4383a7-3ce4-4ecd-8518-c9c7cb681e01",
  "6d9944ec-a347-4b09-8ebc-0e0b2e54875a"
)
experiment_id_name = c (
  "37e41038-310c-44e1-a1ed-8979b0f73709" = "0.000103",
  "082e6a59-27fe-47ed-ae2c-89da8ceb31ee" = "0.000508",
  "ce333fb9-5e8d-4b14-a5ee-7101afe2cfa6" = "0.001001",
  "b3b1b9af-9a6c-43dd-b5ab-be35bbee9c7f" = "0.001017",
  "40a08990-33b3-41f5-b061-2561bc84ca96" = "0.001504",
  "8c689128-8483-4ce5-9b67-8cb160550497" = "0.002613",
  "b7e78bf9-30a3-4955-9469-0bbd77310476" = "0.005028",
  "9442e6ea-8033-44b2-b0f3-91470348a7f0" = "0.006923",
  "156d68c8-abcd-439a-a387-71b3708652ff" = "0.009596",
  "7d4caa7a-5b9c-4b47-af6a-0217ce344f6a" = "0.013481",
  "a362969d-fe30-497a-a2ec-42bc0119f398" = "0.016086",
  "62bdb0f5-affe-4b31-8b27-f5e96899f37f" = "0.019272",
  "b37d56c6-488a-4121-abf6-8b8a2bf82668" = "0.023176",
  "424c4ff7-5cc6-4032-bd54-70b492bbcbb5" = "0.034073",
  "a5b2c471-f1bc-49b7-9650-19bcd78e73ae" = "0.041562",
  "58020739-79fd-4cab-8f4e-e611a6cfba15" = "0.050847",
  "ffad3de5-49d4-4c8d-8be9-6db1be42dbda" = "0.062095",
  "92ffff7e-8c09-447a-96b1-9bcdc2a0da1c" = "0.075168",
  "a11ae384-a142-4f33-a129-eda212675e2a" = "0.108786",
  "44eb6f63-6f19-47d2-acce-67f2e59395bf" = "0.131990",
  "1a0d8c5d-3a1b-48ec-894a-12888c735955" = "0.158308",
  "3c4383a7-3ce4-4ecd-8518-c9c7cb681e01" = "0.191628",
  "6d9944ec-a347-4b09-8ebc-0e0b2e54875a" = "0.231434"
)


configuration_order = c (
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9"
)
configuration_name = c (
  "1" = "deml fastq SI to SS fastq",
  "2" = "pheniqs fastq SI to null",
  "3" = "pheniqs fastq SI to II bam",
  "4" = "pheniqs fastq SI to SI fastq ",
  "5" = "pheniqs fastq SI to SS fastq",
  "6" = "pheniqs fastq SI to II cram",
  "7" = "pheniqs mdd fastq SI to II bam",
  "8" = "bcl2fastq bcl to SI fastq",
  "9" = "bcl2fastq bcl to SS fastq"
)
configuration_color = c (
  "1" = alpha("#C00606", 1),
  "2" = alpha("#5F9A10", 1),
  "3" = alpha("#5F9A10", 1),
  "4" = alpha("#5F9A10", 1),
  "5" = alpha("#5F9A10", 1),
  "6" = alpha("#5F9A10", 1),
  "7" = alpha("#5F9A10", 1),
  "8" = alpha("#3852BB", 1),
  "9" = alpha("#3852BB", 1)
)
configuration_fill_scale = scale_fill_manual (
  name = "Configuration",
  breaks = configuration_order,
  labels = configuration_name,
  values = configuration_color
)

duration_scale = scale_y_continuous (
    breaks = c (
      0,
      # 1800,
      # 3600,
      2 * 3600,
      4 * 3600,
      8 * 3600,
      16 * 3600,
      32 * 3600,
      64 * 3600,
      287509,
      84 * 3600
    ),
    labels = c (
      "00:00:00",
      # "00:30:00",
      # "01:00:00",
      "02:00:00",
      "04:00:00",
      "08:00:00",
      "16:00:00",
      "32:00:00",
      "64:00:00",
      "79:51:49",
      "80:00:00"
    )
)

accurecy_variable_labeller = labeller (
  ssid = experiment_id_name,
  tool = tool_name,
  variable = accurecy_variable_name
)
rate_scale = scale_x_continuous (
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
pheniqs_plot_theme = theme (
  plot.title = plot_title_text,
  plot.margin = margin(0.5,0.5,0.5,0.5, "char"),
  axis.line = axis_line,
  axis.ticks = tick_line,
  axis.title.x = title_text,
  axis.text.x = vertical_axis_text,
  axis.title.y = vertical_axis_text,
  axis.text.y = axis_text,
  strip.text = title_text,
  strip.background = element_blank(),
  strip.placement = "outside",
  panel.grid.major = grid_line,
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = legend_title_text,
  legend.text = legend_text,
  legend.key.width = unit(0.375, "char"),
  legend.key.height = unit(0.375, "char"),
  legend.key = element_rect(fill = NA),
  legend.margin = margin(0,0,0,0, "char")
)

draw_diagram <- function(plot, diagram_filename, diagram_width, diagram_height) {
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
}
