#!/usr/bin/env zsh

PHENIQS_HOME=~/code/moonwatcher/pheniqs
REMOTE_SIMULATION_HOME=albireo:/Volumes/canal/A5KVK_3
REMOTE_HOME=Sites/secret/2d1412bc-9643-44de-889e-346a427370a4/A5KVK_V3

MANUSCRIPT_HOME=$PHENIQS_HOME/manuscript
PLOTTING_CODE_DIRECTORY=$MANUSCRIPT_HOME/plotting
MANUSCRIPT_CSV_DIRECTORY=$MANUSCRIPT_HOME/csv
MANUSCRIPT_PLOT_DIRECTORY=$MANUSCRIPT_HOME/plot
SIMULATION_HOME=$MANUSCRIPT_HOME

BARCODE_SIMULATION_ID=2ea833ea-d9c5-4994-a4f3-d4f786e7e19a
SIMULATION_SESSION_ID=2ea833ea

clean() {
    echo "Cleaning..."
    if [ -d $MANUSCRIPT_CSV_DIRECTORY ]; then
        echo "Cleaning csv"
        rm -rf $MANUSCRIPT_CSV_DIRECTORY
    fi

    if [ -d $MANUSCRIPT_PLOT_DIRECTORY ]; then
        echo "Cleaning plot"
        rm -rf $MANUSCRIPT_PLOT_DIRECTORY
    fi
}

prepare() {
    if [ ! -d $SIMULATION_HOME ]; then
        echo "Preparing simulation home directory"
        mkdir -p $SIMULATION_HOME
    fi

    if [ ! -d $MANUSCRIPT_CSV_DIRECTORY ]; then
        echo "Preparing csv directory"
        mkdir -p $MANUSCRIPT_CSV_DIRECTORY
    fi

    if [ ! -d $MANUSCRIPT_PLOT_DIRECTORY ]; then
        echo "Preparing plot directory"
        mkdir -p $MANUSCRIPT_PLOT_DIRECTORY
    fi
}

fetch_remote() {
    rsync --progress $REMOTE_SIMULATION_HOME/benchmark.json.bz2 $SIMULATION_HOME/
    rsync --recursive --progress --exclude '*.bam' $REMOTE_SIMULATION_HOME/simulation $SIMULATION_HOME/
}

extract_archive() {
    if [ ! -f $SIMULATION_HOME/benchmark.json ]; then
        echo "Extracting simulation archive"
        tar -xvjf simulation.tar.bz2
    fi
}

make_csv() {
  PRESET=$1
  FILENAME=$2.csv

  if [ ! -f $MANUSCRIPT_CSV_DIRECTORY/$FILENAME ]; then
      echo "Compiling $FILENAME"
      $PHENIQS_HOME/tool/pheniqs-benchmark-api.py \
      --home $SIMULATION_HOME \
      summarize \
      --preset $PRESET \
      --bsid $BARCODE_SIMULATION_ID \
      --session $SIMULATION_SESSION_ID \
      > $MANUSCRIPT_CSV_DIRECTORY/$FILENAME;
  fi
};

make_plot() {
  CSV_FILENAME=$1.csv
  R_FILENAME=$2.R
  PLOT_FILENAME=$2.pdf

  if [ ! -f $MANUSCRIPT_PLOT_DIRECTORY/$PLOT_FILENAME ]; then
      echo "Plotting $PLOT_FILENAME"
      ( cd $PLOTTING_CODE_DIRECTORY;
        ./$R_FILENAME \
        $MANUSCRIPT_CSV_DIRECTORY/$CSV_FILENAME \
        $MANUSCRIPT_PLOT_DIRECTORY/$PLOT_FILENAME;

        [ -f Rplots.pdf ] && rm -f Rplots.pdf
      )
  fi
};

make_data() {
    make_csv decoder_summary_R decoder_summary_R
    make_csv decoder_summary decoder_summary
    make_csv barcode_summary_R barcode_summary_R
    make_csv barcode_summary barcode_summary
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
    make_csv binned_barcode_prior binned_barcode_prior
};

make_plots() {
    # 0550 prefix
    # maximum expected substitution error 0.055
    # 5.5% substituted nucleotides
    # 55 substituted nucleotides in a 1000

    # 0060 prefix
    # maximum expected substitution error 0.006
    # 0.6% substituted nucleotides
    # 6 substituted nucleotides in a 1000

    make_plot barcode_distribution barcode_distribution
    make_plot classifiable_summary_R classifiable_accuracy_0550
    make_plot classified_summary_R classified_accuracy_0060_mdd
    make_plot binned_decoder_summary_R classified_accuracy_0550_binned
    make_plot classified_summary_R classified_accuracy_0550
    make_plot decoder_summary_R overall_0550_mdd
    make_plot decoder_summary_R overall_0550
    make_plot binned_barcode_prior prior_estimation_error_0550_binned
    make_plot quality_distribution quality_distribution_by_rate_0550
    make_plot noise_summary_R unclassified_accuracy_0060_mdd
    make_plot noise_summary_R unclassified_accuracy_0550
};

make_all_plot() {
    make_plot barcode_distribution barcode_prior_distribution
    make_plot quality_distribution quality_distribution_by_rate_0550
    make_plot binned_barcode_prior prior_estimation_error_binned_0550
    make_plot classifiable_summary_R classifiable_accuracy_0550
    # make_plot classifiable_summary_R classifiable_overview_0550_mdd
    make_plot classified_summary_R classified_accuracy_0550
    # make_plot classified_summary_R classified_overview_0550_mdd
    make_plot decoder_summary_R overall_0550
    make_plot decoder_summary_R overall_0550_mdd
    make_plot noise_summary_R unclassified_accuracy_0550
    # make_plot noise_summary_R unclassified_overview_0550_mdd
    make_plot binned_decoder_summary_R classified_binned_0550
    # make_plot binned_decoder_summary_R binned_classified_0550_mdd



    # maximum error rate 0.055 / 5.5% / 55 in a 1000

    make_plot barcode_distribution barcode_true_prior
    make_plot quality_distribution quality_distribution_0550
    make_plot binned_barcode_prior binned_prior_0550
    make_plot classifiable_summary_R classifiable_overview_0550
    # make_plot classifiable_summary_R classifiable_overview_0550_mdd
    make_plot classified_summary_R classified_overview_0550
    # make_plot classified_summary_R classified_overview_0550_mdd
    make_plot decoder_summary_R overall_overview_0550
    make_plot decoder_summary_R overall_overview_0550_mdd
    make_plot noise_summary_R unclassified_overview_0550
    # make_plot noise_summary_R unclassified_overview_0550_mdd
    make_plot binned_decoder_summary_R binned_classified_0550
    # make_plot binned_decoder_summary_R binned_classified_0550_mdd

    # maximum error rate 0.006 / 0.6% / 6 in a 1000
    #
    # make_plot quality_distribution quality_distribution_0060
    # make_plot binned_barcode_prior binned_prior_0060
    # make_plot classifiable_summary_R classifiable_overview_0060
    # make_plot classifiable_summary_R classifiable_overview_0060_mdd
    # make_plot classified_summary_R classified_overview_0060
    make_plot classified_summary_R classified_overview_0060_mdd
    # make_plot decoder_summary_R overall_overview_0060
    # make_plot decoder_summary_R overall_overview_0060_mdd
    # make_plot noise_summary_R unclassified_overview_0060
    make_plot noise_summary_R unclassified_overview_0060_mdd
    # make_plot binned_decoder_summary_R binned_classified_0060
    # make_plot binned_decoder_summary_R binned_classified_0060_mdd

    # make_plot quality_distribution quality_distribution_0027
    # make_plot binned_barcode_prior binned_prior_0027
    # make_plot classifiable_summary_R classifiable_overview_0027
    # make_plot classifiable_summary_R classifiable_overview_0027_mdd
    # make_plot classified_summary_R classified_overview_0027
    # make_plot classified_summary_R classified_overview_0027_mdd
    # make_plot decoder_summary_R overall_overview_0027
    # make_plot decoder_summary_R overall_overview_0027_mdd
    # make_plot noise_summary_R unclassified_overview_0027
    # make_plot noise_summary_R unclassified_overview_0027_mdd
    # make_plot binned_decoder_summary_R binned_classified_0027
    # make_plot binned_decoder_summary_R binned_classified_0027_mdd

}

push_remote () {
    rsync --recursive --progress $MANUSCRIPT_CSV_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
    rsync --recursive --progress $MANUSCRIPT_PLOT_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
}

usage() {
    echo '-f --fetch     fetch remote simulation'
    echo '-p --push      push remote simulation'
    echo '-c --clean     clean csv and plot directory'
    echo '-h --help      print this help'
}

clean_set=0
fetch_remote_set=0
push_remote_set=0
while [ "$1" != "" ]; do
    case $1 in
        -c | --clean )          clean_set=1
                                ;;
        -f | --fetch )          fetch_remote_set=1
                                ;;
        -p | --push )           push_remote_set=1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

[ ! $clean_set = 0 ] &&  clean;

prepare

[ ! $fetch_remote_set = 0 ] && fetch_remote;

extract_archive
make_data
make_plots

[ ! $push_remote_set = 0 ] && push_remote;
