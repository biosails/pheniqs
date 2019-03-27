#!/usr/bin/env zsh

PHENIQS_HOME=~/code/biosails/pheniqs

EXPERIMENT_HOME=$PHENIQS_HOME/local/A5KVK
[ -d $EXPERIMENT_HOME ] || mkdir -p $EXPERIMENT_HOME

# [ -d $EXPERIMENT_HOME/data ] && rm -rf $EXPERIMENT_HOME/data
[ -d $EXPERIMENT_HOME/data ] || mkdir -p $EXPERIMENT_HOME/data
# benchmark.py summarize --session 2ea833ea -p decoder_summary_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/decoder_summary_R.csv
# benchmark.py summarize --session 2ea833ea -p decoder_summary --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/decoder_summary.csv
# benchmark.py summarize --session 2ea833ea -p barcode_summary_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/barcode_summary_R.csv
# benchmark.py summarize --session 2ea833ea -p barcode_summary --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/barcode_summary.csv
# benchmark.py summarize --session 2ea833ea -p noise_summary_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/noise_summary_R.csv
# benchmark.py summarize --session 2ea833ea -p noise_summary --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/noise_summary.csv
# benchmark.py summarize --session 2ea833ea -p classified_summary_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/classified_summary_R.csv
# benchmark.py summarize --session 2ea833ea -p classified_summary --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/classified_summary.csv
# benchmark.py summarize --session 2ea833ea -p binned_decoder_summary_R --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/binned_decoder_summary_R.csv
# benchmark.py summarize --session 2ea833ea -p binned_decoder_summary --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/binned_decoder_summary.csv
# benchmark.py summarize --session 2ea833ea -p quality_distribution --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/quality_distribution.csv
# benchmark.py summarize --session 2ea833ea -p barcode_distribution --bsid 2ea833ea-d9c5-4994-a4f3-d4f786e7e19a > $EXPERIMENT_HOME/data/barcode_distribution.csv
#
[ -d $EXPERIMENT_HOME/plot ] && rm -rf $EXPERIMENT_HOME/plot
[ -d $EXPERIMENT_HOME/plot ] || mkdir -p $EXPERIMENT_HOME/plot
( cd $PHENIQS_HOME/tool/simulation/plot;

./plot_decoder_accuracy.R $EXPERIMENT_HOME/data/decoder_summary_R.csv $EXPERIMENT_HOME/plot/decoder_summary.pdf
./plot_noise_accuracy.R $EXPERIMENT_HOME/data/noise_summary_R.csv $EXPERIMENT_HOME/plot/noise_summary.pdf
./plot_classified_accuracy.R $EXPERIMENT_HOME/data/classified_summary_R.csv $EXPERIMENT_HOME/plot/classified_summary.pdf
./plot_binned_accuracy.R $EXPERIMENT_HOME/data/binned_decoder_summary_R.csv $EXPERIMENT_HOME/plot/binned_decoder_summary.pdf
./plot_comparison.R $EXPERIMENT_HOME/data/classified_summary.csv $EXPERIMENT_HOME/plot/comparison.pdf
./plot_quality_distribution.R $EXPERIMENT_HOME/data/quality_distribution.csv $EXPERIMENT_HOME/plot/quality_distribution.pdf
./plot_barcode_distribution.R $EXPERIMENT_HOME/data/barcode_distribution.csv $EXPERIMENT_HOME/plot/barcode_distribution.pdf
./plot_fscore.R $EXPERIMENT_HOME/data/decoder_summary_R.csv $EXPERIMENT_HOME/plot/decoder_fscore.pdf

rm Rplots.pdf )

rsync --recursive --progress $EXPERIMENT_HOME/data lg@albireo.bio.nyu.edu:Sites/secret/2d1412bc-9643-44de-889e-346a427370a4/A5KVK/
rsync --recursive --progress $EXPERIMENT_HOME/plot lg@albireo.bio.nyu.edu:Sites/secret/2d1412bc-9643-44de-889e-346a427370a4/A5KVK/
