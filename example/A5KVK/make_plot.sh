#!/usr/bin/env zsh


PHENIQS_HOME=~/code/biosails/pheniqs

REMOTE_SIMULATION_HOME=albireo:/Volumes/canal/A5KVK_3
SIMULATION_HOME=$PHENIQS_HOME/example/A5KVK
SIMULATION_CSV_DIRECTORY=$SIMULATION_HOME/csv
SIMULATION_PLOT_DIRECTORY=$SIMULATION_HOME/plot
BARCODE_SIMULATION_ID=2ea833ea-d9c5-4994-a4f3-d4f786e7e19a
SIMULATION_SESSION_ID=2ea833ea
# MAXIMUM_ERROR_RATE=0.076
REMOTE_HOME=Sites/secret/2d1412bc-9643-44de-889e-346a427370a4/A5KVK_V3

prepare() {
    if [ ! -d $SIMULATION_CSV_DIRECTORY ]; then
        echo "Preparing csv directory"
        mkdir -p $SIMULATION_CSV_DIRECTORY
    fi

    if [ ! -d $SIMULATION_PLOT_DIRECTORY ]; then
        echo "Preparing plot directory"
        mkdir -p $SIMULATION_PLOT_DIRECTORY
    fi
}

clean() {
    if [ -d $SIMULATION_CSV_DIRECTORY ]; then
        echo "Cleaning csv"
        rm -rf $SIMULATION_CSV_DIRECTORY
    fi

    if [ -d $SIMULATION_PLOT_DIRECTORY ]; then
        echo "Cleaning plot"
        rm -rf $SIMULATION_PLOT_DIRECTORY
    fi
}

extract_archive() {
    if [ ! -f $SIMULATION_HOME/benchmark.json ]; then
        echo "Extracting benchmark.json"
        bzip2 -d -k benchmark.json.bz2
    fi
}

make_csv() {
  PRESET=$1
  CSV_PATH=$2.csv

  if [ ! -f $SIMULATION_CSV_DIRECTORY/$CSV_PATH ]; then
      $PHENIQS_HOME/tool/pheniqs-benchmark-api.py \
      --home $SIMULATION_HOME \
      summarize \
      --bsid $BARCODE_SIMULATION_ID \
      --session $SIMULATION_SESSION_ID \
      --preset $PRESET \
      > $SIMULATION_CSV_DIRECTORY/$CSV_PATH;
  fi
};

make_plot() {
  R_METHOD_PATH=$1.R
  CSV_PATH=$2.csv
  MAXIMUM_ERROR_RATE=0.$4
  PLOT_PATH=$3_$5$4.pdf

  if [ ! -f $SIMULATION_PLOT_DIRECTORY/$PLOT_PATH ]; then
      ( cd $PHENIQS_HOME/tool/simulation/plot;
        ./$R_METHOD_PATH \
        $SIMULATION_CSV_DIRECTORY/$CSV_PATH \
        $SIMULATION_PLOT_DIRECTORY/$PLOT_PATH \
        $MAXIMUM_ERROR_RATE \
        $5;

        [ -f Rplots.pdf ] && rm -f Rplots.pdf
      )
  fi
};

make_all_csv() {
    make_csv decoder_summary_R decoder_summary_R
    make_csv decoder_summary decoder_summary
    make_csv barcode_summary_R barcode_summary_R
    make_csv barcode_summary barcode_summary
    make_csv barcode_prior barcode_prior
    make_csv noise_summary_R noise_summary_R
    make_csv noise_summary noise_summary
    make_csv classified_summary_R classified_summary_R
    make_csv classified_summary classified_summary
    make_csv classifiable_summary_R classifiable_summary_R
    make_csv classifiable_summary classifiable_summary
    make_csv binned_decoder_summary_R binned_decoder_summary_R
    make_csv binned_decoder_summary binned_decoder_summary
    make_csv quality_distribution quality_distribution
    make_csv barcode_distribution barcode_distribution
}

make_all_plot() {
    ERROR_LIMIT=$1
    no_mdd=$2
    make_plot plot_decoder_accuracy decoder_summary_R decoder_summary $ERROR_LIMIT $no_mdd
    make_plot plot_noise_accuracy noise_summary_R noise_summary $ERROR_LIMIT $no_mdd
    make_plot plot_classified_accuracy classified_summary_R classified_summary $ERROR_LIMIT $no_mdd
    make_plot plot_classifiable_accuracy classifiable_summary_R classifiable_summary $ERROR_LIMIT $no_mdd
    make_plot plot_binned_accuracy binned_decoder_summary_R binned_decoder_summary $ERROR_LIMIT $no_mdd
    make_plot plot_comparison classified_summary comparison $ERROR_LIMIT $no_mdd
    make_plot plot_quality_distribution quality_distribution quality_distribution $ERROR_LIMIT $no_mdd
    make_plot plot_barcode_distribution barcode_distribution barcode_distribution $ERROR_LIMIT $no_mdd
    make_plot plot_decoder_prior barcode_prior decoder_prior $ERROR_LIMIT $no_mdd
    make_plot plot_barcode_prior barcode_prior barcode_prior $ERROR_LIMIT $no_mdd

    make_plot plot_fscore decoder_summary_R decoder_fscore $ERROR_LIMIT $no_mdd
    make_plot plot_precision_recall decoder_summary_R precision_recall $ERROR_LIMIT $no_mdd
    # make_plot plot_FDR_MR decoder_summary_R FDR_MR $ERROR_LIMIT $no_mdd
    make_plot plot_classified_FDR_MR classified_summary_R classified_FDR_MR $ERROR_LIMIT $no_mdd
    make_plot plot_classifiable_FDR_MR classifiable_summary_R classifiable_FDR_MR $ERROR_LIMIT $no_mdd
}

fetch_remote() {
    rsync --progress $REMOTE_SIMULATION_HOME/benchmark.json.bz2 $SIMULATION_HOME/
    rsync --recursive --progress --exclude '*.bam' $REMOTE_SIMULATION_HOME/simulation $SIMULATION_HOME/
}

push_remote () {
    rsync --recursive --progress $SIMULATION_CSV_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
    rsync --recursive --progress $SIMULATION_PLOT_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
}

# verify the simulation home directory exists
[ -d $SIMULATION_HOME ] || mkdir -p $SIMULATION_HOME

# clean output directories for csv and plot
# clean

# verify output directories for csv and plot files exist
prepare

# fetch simulation data from remote location
# fetch_remote

# extract the archive if not present
extract_archive

# generate csv files
make_all_csv

# generate plot files
make_all_plot 0027
make_all_plot 0027 no_mdd
make_all_plot 01
make_all_plot 01 no_mdd
make_all_plot 03
make_all_plot 03 no_mdd
make_all_plot 076
make_all_plot 076 no_mdd
make_all_plot 2
make_all_plot 2 no_mdd

# push csv and plot to remote location
# push_remote
