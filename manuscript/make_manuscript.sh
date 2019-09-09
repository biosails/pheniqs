#!/usr/bin/env zsh

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2018  Lior Galanti
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

PHENIQS_HOME=~/code/moonwatcher/pheniqs
REMOTE_SIMULATION_HOME=albireo:/Volumes/canal/A5KVK_3
REMOTE_HOME=Sites/secret/2d1412bc-9643-44de-889e-346a427370a4/A5KVK_V3

MANUSCRIPT_HOME=$PHENIQS_HOME/manuscript
PLOTTING_CODE_DIRECTORY=$MANUSCRIPT_HOME/plot_script
MANUSCRIPT_CSV_DIRECTORY=$MANUSCRIPT_HOME/csv
MANUSCRIPT_PLOT_DIRECTORY=$MANUSCRIPT_HOME/plot
MANUSCRIPT_BIORXIV_DIRECTORY=$MANUSCRIPT_HOME/biorxiv
MANUSCRIPT_BMC_DIRECTORY=$MANUSCRIPT_HOME/bmc
SIMULATION_HOME=$MANUSCRIPT_HOME

BARCODE_SIMULATION_ID=2ea833ea-d9c5-4994-a4f3-d4f786e7e19a
SIMULATION_SESSION_ID=2ea833ea

MANUSCRIPT_NAME="pheniqs"

make_simulation_csv() {
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

make_r_plot() {
  CSV_FILENAME=$1.csv
  R_FILENAME=$2.R
  PLOT_FILENAME=$2.pdf
  [ $3 ] && PLOT_FILENAME=$3_$PLOT_FILENAME

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

push_remote () {
    rsync --recursive --progress $MANUSCRIPT_CSV_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
    rsync --recursive --progress $MANUSCRIPT_PLOT_DIRECTORY lg@albireo.bio.nyu.edu:$REMOTE_HOME/
}

fetch_remote() {
    rsync --progress $REMOTE_SIMULATION_HOME/benchmark.json.bz2 $SIMULATION_HOME/
    rsync --recursive --progress --exclude '*.bam' $REMOTE_SIMULATION_HOME/simulation $SIMULATION_HOME/
}

clean_all() {
    if [ -d $MANUSCRIPT_CSV_DIRECTORY ]; then
        echo "Removing csv directory"
        rm -rf $MANUSCRIPT_CSV_DIRECTORY
    fi

    if [ -d $MANUSCRIPT_PLOT_DIRECTORY ]; then
        echo "Removing plot directory"
        rm -rf $MANUSCRIPT_PLOT_DIRECTORY
    fi

    if [ -d $MANUSCRIPT_HOME/simulation ]; then
        echo "Removing simulation directory"
        rm -rf $MANUSCRIPT_HOME/simulation
    fi

    if [ -d $MANUSCRIPT_HOME/session ]; then
        echo "Removing session directory"
        rm -rf $MANUSCRIPT_HOME/session
    fi

    if [ -f $MANUSCRIPT_HOME/benchmark.json ]; then
        echo "Removing benchmark.json"
        rm -rf $MANUSCRIPT_HOME/benchmark.json
    fi

    if [ -f $MANUSCRIPT_HOME/benchmark.pickle ]; then
        echo "Removing benchmark.pickle"
        rm -rf $MANUSCRIPT_HOME/benchmark.pickle
    fi
}

prepare() {
    if [ -d $MANUSCRIPT_CSV_DIRECTORY ]; then
        if [ ! $clean_csv_set = 0 ]; then
            echo "Removing existing csv directory"
            rm -rf $MANUSCRIPT_CSV_DIRECTORY
            clean_csv_set=0
            clean_plot_set=1
            make_csv_set=1
            make_plot_set=1
        fi
    fi

    # clean existing directories
    if [ -d $MANUSCRIPT_PLOT_DIRECTORY ]; then
        if [ ! $clean_plot_set = 0 ]; then
            echo "Removing existing plot directory"
            rm -rf $MANUSCRIPT_PLOT_DIRECTORY
            make_plot_set=1
        fi
    fi

    if [ -f $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_biorxiv.pdf ]; then
        if [ ! $clean_biorxiv_set = 0 ]; then
            echo "Removing biorxiv manuscript"
            rm $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_biorxiv.pdf
            make_biorxiv_set=1
        fi
    fi

    if [ -f $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_bmc.pdf ]; then
        if [ ! $clean_bmc_set = 0 ]; then
            echo "Removing bmc manuscript"
            rm $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_bmc.pdf
            make_bmc_set=1
        fi
    fi

    # varifying directories
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
    if [ ! -d $MANUSCRIPT_BIORXIV_DIRECTORY ]; then
        echo "Preparing biorxiv directory"
        mkdir -p $MANUSCRIPT_BIORXIV_DIRECTORY
    fi
    if [ ! -d $MANUSCRIPT_BMC_DIRECTORY ]; then
        echo "Preparing bmc directory"
        mkdir -p $MANUSCRIPT_BMC_DIRECTORY
    fi

    # expand archive
    if [ ! -f $SIMULATION_HOME/benchmark.json ]; then
        echo "Extracting simulation archive"
        tar -xjf simulation.tar.bz2
    fi
}

make_csv() {
    make_simulation_csv decoder_summary_R decoder_summary_R
    make_simulation_csv decoder_summary decoder_summary
    make_simulation_csv barcode_summary_R barcode_summary_R
    make_simulation_csv barcode_summary barcode_summary
    make_simulation_csv noise_summary_R noise_summary_R
    make_simulation_csv noise_summary noise_summary
    make_simulation_csv classified_summary_R classified_summary_R
    make_simulation_csv classified_summary classified_summary
    make_simulation_csv classifiable_summary_R classifiable_summary_R
    make_simulation_csv classifiable_summary classifiable_summary
    make_simulation_csv binned_decoder_summary_R binned_decoder_summary_R
    make_simulation_csv binned_decoder_summary binned_decoder_summary
    make_simulation_csv quality_distribution quality_distribution
    make_simulation_csv barcode_distribution barcode_distribution
    make_simulation_csv binned_barcode_prior binned_barcode_prior
    make_csv_set=0
};

make_plot() {
    make_csv
    # 0550 prefix
    # maximum expected substitution error 0.055
    # 5.5% substituted nucleotides
    # 55 substituted nucleotides in a 1000

    # 0060 prefix
    # maximum expected substitution error 0.006
    # 0.6% substituted nucleotides
    # 6 substituted nucleotides in a 1000

    make_r_plot decoder_summary_R           overall_0550_mdd                    1
    make_r_plot decoder_summary_R           overall_0550                        2
    make_r_plot classifiable_summary_R      classifiable_accuracy_0550          3
    make_r_plot classifiable_summary_R      classifiable_accuracy_0060          4
    make_r_plot classified_summary_R        classified_accuracy_0550            5
    make_r_plot classified_summary_R        classified_accuracy_0060_mdd        6
    make_r_plot binned_decoder_summary_R    classified_accuracy_0550_binned     7
    make_r_plot binned_decoder_summary_R    classified_accuracy_0060_binned     8
    make_r_plot noise_summary_R             unclassified_accuracy_0550          9
    make_r_plot noise_summary_R             unclassified_accuracy_0060_mdd      10
    make_r_plot noise_summary_R             unclassified_accuracy_0060          11
    make_r_plot binned_barcode_prior        prior_estimation_error_0550_binned  12
    make_r_plot barcode_distribution        barcode_distribution                13
    make_r_plot quality_distribution        quality_distribution_by_rate_0550   14
    make_r_plot ../speed_memory             speed_memory                        15
    make_plot_set=0
};

make_biorxiv() {
    if [ ! -f $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_biorxiv.pdf ]; then
        make_plot
        [ -f $MANUSCRIPT_BIORXIV_DIRECTORY/$MANUSCRIPT_NAME.pdf ] && rm $MANUSCRIPT_BIORXIV_DIRECTORY/$MANUSCRIPT_NAME.pdf

        (   cd "$MANUSCRIPT_BIORXIV_DIRECTORY";
            [[ -f $MANUSCRIPT_NAME.aux ]] && rm $MANUSCRIPT_NAME.aux;
            [[ -f $MANUSCRIPT_NAME.log ]] && rm $MANUSCRIPT_NAME.log;
            [[ -f $MANUSCRIPT_NAME.bbl ]] && rm $MANUSCRIPT_NAME.bbl;
            [[ -f $MANUSCRIPT_NAME.blg ]] && rm $MANUSCRIPT_NAME.blg;
            [[ -f $MANUSCRIPT_NAME.bcf ]] && rm $MANUSCRIPT_NAME.bcf;
            [[ -f $MANUSCRIPT_NAME.run.xml ]] && rm $MANUSCRIPT_NAME.run.xml;
            [[ -f $MANUSCRIPT_NAME.pdf ]] && rm $MANUSCRIPT_NAME.pdf;

            pdflatex $MANUSCRIPT_NAME > /dev/null && \
            biber $MANUSCRIPT_NAME > /dev/null && \
            pdflatex $MANUSCRIPT_NAME > /dev/null && \
            pdflatex $MANUSCRIPT_NAME > /dev/null;

            [[ -f $MANUSCRIPT_NAME.aux ]] && rm $MANUSCRIPT_NAME.aux;
            [[ -f $MANUSCRIPT_NAME.log ]] && rm $MANUSCRIPT_NAME.log;
            [[ -f $MANUSCRIPT_NAME.bbl ]] && rm $MANUSCRIPT_NAME.bbl;
            [[ -f $MANUSCRIPT_NAME.blg ]] && rm $MANUSCRIPT_NAME.blg;
            [[ -f $MANUSCRIPT_NAME.bcf ]] && rm $MANUSCRIPT_NAME.bcf;
            [[ -f $MANUSCRIPT_NAME.run.xml ]] && rm $MANUSCRIPT_NAME.run.xml;
        )
        cp $MANUSCRIPT_BIORXIV_DIRECTORY/$MANUSCRIPT_NAME.pdf $MANUSCRIPT_HOME/${MANUSCRIPT_NAME}_biorxiv.pdf
        make_biorxiv_set=0
    fi
}

make_bmc() {

    make_bmc_set=0
}

usage() {
    echo '-p --plot             build plot directory'
    echo '-P --clean-plot       rebuild plot directory'
    echo '-c --csv              build csv directory'
    echo '-C --clean-csv        rebuild csv directory'
    echo '-x --biorxiv          build biorxiv manuscript'
    echo '-X --rebuild-biorxiv  rebuild biorxiv manuscript'
    echo '-b --bmc              build bmc manuscript'
    echo '-B --rebuild-bmc      rebuild bmc manuscript'
    echo '-e --clean            remove all generated files'
    # echo '-f --fetch      fetch remote simulation'
    # echo '-p --push       push remote simulation'
    echo '-h --help           print this help'
}

make_plot_set=0
clean_plot_set=0
make_csv_set=0
clean_csv_set=0
make_biorxiv_set=0
clean_biorxiv_set=0
make_bmc_set=0
clean_bmc_set=0
clean_all_set=0

while [ "$1" != "" ]; do
    case $1 in
        -p | --plot )           make_plot_set=1
                                ;;
        -P | --rebuild-plot )   clean_plot_set=1
                                ;;
        -c | --csv )            make_csv_set=1
                                ;;
        -C | --rebuild-csv )    clean_csv_set=1
                                ;;
        -x | --biorxiv )        make_biorxiv_set=1
                                ;;
        -X | --rebuild-biorxiv ) clean_biorxiv_set=1
                                ;;
        -c | --bmc )            make_bmc_set=1
                                ;;
        -C | --rebuild-bmc )    clean_bmc_set=1
                                ;;
        -e | --clean )          clean_all_set=1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

if [ ! $clean_all_set = 0 ]; then
    clean_all
else
    prepare

    if [ ! $fetch_remote_set = 0 ]; then
        fetch_remote
    fi

    if [ ! $make_plot_set = 0 ]; then
        [ ! $make_csv_set = 0 ] && make_csv
        make_plot
    fi

    if [ ! $make_csv_set = 0 ]; then
        make_csv
    fi

    if [ ! $make_biorxiv_set = 0 ]; then
        [ ! $make_plot_set = 0 ] && make_plot
        make_biorxiv
    fi

    if [ ! $make_bmc_set = 0 ]; then
        [ ! $make_plot_set = 0 ] && make_plot
        make_bmc
    fi

    if [ ! $push_remote_set = 0 ]; then
        push_remote
    fi
fi
