#!/usr/bin/env zsh

PHENIQS_HOME=~/code/biosails/pheniqs

# small original
# benchmark.py -H ~/Downloads/A5KVK.original summarize -p multiplex --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex.csv
# benchmark.py -H ~/Downloads/A5KVK.original summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex_R.csv
# benchmark.py -H ~/Downloads/A5KVK.original summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > noise_R.csv
# benchmark.py -H ~/Downloads/A5KVK.original summarize -p quality_distribution_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > quality_distribution_R.csv

# small alternative
# benchmark.py -H ~/Downloads/A5KVK.alternative summarize -p multiplex --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex.csv
# benchmark.py -H ~/Downloads/A5KVK.alternative summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex_R.csv
# benchmark.py -H ~/Downloads/A5KVK.alternative summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > noise_R.csv
# benchmark.py -H ~/Downloads/A5KVK.alternative summarize -p quality_distribution_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > quality_distribution_R.csv

# small
# benchmark.py -H ~/Downloads/A5KVK summarize -p multiplex --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex.csv
# benchmark.py -H ~/Downloads/A5KVK summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex_R.csv
# benchmark.py -H ~/Downloads/A5KVK summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > noise_R.csv
# benchmark.py -H ~/Downloads/A5KVK summarize -p quality_distribution_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > quality_distribution_R.csv

# second
# benchmark.py summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > multiplex_R.csv
# benchmark.py summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > noise_R.csv

# first
# EXPERIMENT_HOME=$PHENIQS_HOME/local/A5KVK/first/plot
# [ -d $EXPERIMENT_HOME ] || mkdir -p $EXPERIMENT_HOME
#
# benchmark.py summarize -p multiplex --bsid 3745482b-c557-42bb-a842-565e3b59a308 > $EXPERIMENT_HOME/multiplex.csv
# benchmark.py summarize -p multiplex_R --bsid 3745482b-c557-42bb-a842-565e3b59a308 > $EXPERIMENT_HOME/multiplex_R.csv
# benchmark.py summarize -p noise_R --bsid 3745482b-c557-42bb-a842-565e3b59a308 > $EXPERIMENT_HOME/noise_R.csv
# benchmark.py summarize -p quality_distribution_R --bsid 3745482b-c557-42bb-a842-565e3b59a308 > $EXPERIMENT_HOME/quality_distribution_R.csv
#
# $PHENIQS_HOME/tool/simulation/plot/plot_real_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/real_multiplex.pdf
# $PHENIQS_HOME/tool/simulation/plot/plot_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/multiplex.pdf
# $PHENIQS_HOME/tool/simulation/plot/plot_noise_accuracy.R $EXPERIMENT_HOME/noise_R.csv $EXPERIMENT_HOME/noise.pdf
# $PHENIQS_HOME/tool/simulation/plot/plot_quality_distribution.R $EXPERIMENT_HOME/quality_distribution_R.csv $EXPERIMENT_HOME/quality_distribution.pdf
# $PHENIQS_HOME/tool/simulation/plot/plot_comparison.R $EXPERIMENT_HOME/multiplex.csv $EXPERIMENT_HOME/comparison.pdf

# second
EXPERIMENT_HOME=$PHENIQS_HOME/local/A5KVK/plot
[ -d $EXPERIMENT_HOME ] || mkdir -p $EXPERIMENT_HOME

rm $EXPERIMENT_HOME/*.pdf
rm $EXPERIMENT_HOME/*.csv

benchmark.py -v debug summarize -p multiplex --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/multiplex.csv
benchmark.py -v debug summarize -p multiplex_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/multiplex_R.csv
benchmark.py -v debug summarize -p noise_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/noise_R.csv
benchmark.py -v debug summarize -p quality_distribution_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/quality_distribution_R.csv

$PHENIQS_HOME/tool/simulation/plot/plot_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/multiplex.pdf
$PHENIQS_HOME/tool/simulation/plot/plot_real_accuracy.R $EXPERIMENT_HOME/multiplex_R.csv $EXPERIMENT_HOME/real_multiplex.pdf
$PHENIQS_HOME/tool/simulation/plot/plot_noise_accuracy.R $EXPERIMENT_HOME/noise_R.csv $EXPERIMENT_HOME/noise.pdf
$PHENIQS_HOME/tool/simulation/plot/plot_quality_distribution.R $EXPERIMENT_HOME/quality_distribution_R.csv $EXPERIMENT_HOME/quality_distribution.pdf
$PHENIQS_HOME/tool/simulation/plot/plot_comparison.R $EXPERIMENT_HOME/multiplex.csv $EXPERIMENT_HOME/comparison.pdf

 rm Rplots.pdf
