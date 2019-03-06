#!/usr/bin/env zsh

PHENIQS_HOME=~/code/biosails/pheniqs

EXPERIMENT_HOME=$PHENIQS_HOME/local/A5KVK/plot
[ -d $EXPERIMENT_HOME ] || mkdir -p $EXPERIMENT_HOME

rm $EXPERIMENT_HOME/*.pdf
rm $EXPERIMENT_HOME/*.csv

benchmark.py summarize -p multiplex --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/multiplex.csv
benchmark.py summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/multiplex_R.csv
benchmark.py summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/noise_R.csv
benchmark.py summarize -p quality_distribution_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/quality_distribution_R.csv

(
  cd $PHENIQS_HOME/tool/simulation/plot;
  ./plot_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/multiplex.pdf
  ./plot_real_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/real_multiplex.pdf
  ./plot_noise_accuracy.R $EXPERIMENT_HOME/noise_R.csv $EXPERIMENT_HOME/noise.pdf
  ./plot_quality_distribution.R $EXPERIMENT_HOME/quality_distribution_R.csv $EXPERIMENT_HOME/quality_distribution.pdf
  ./plot_comparison.R $EXPERIMENT_HOME/multiplex.csv $EXPERIMENT_HOME/comparison.pdf
  rm Rplots.pdf
)
