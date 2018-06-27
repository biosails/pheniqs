#!/bin/zsh

PROJECT_HOME="$HOME/code/pheniqs";
SIMULATION_DIR="$PROJECT_HOME/benchmark/HK5NHBGXX/simulation"
BENCHMAKR_DIR="$PROJECT_HOME/benchmark";

benchmark_plot() {
    PDF_NAME="benchmark.pdf";
    (   cd "$BENCHMAKR_DIR";
        [[ -f "$PDF_NAME" ]] && rm "$PDF_NAME"
        cat benchmark.log | ./benchmark.py > benchmark.csv;

        ./benchmark.R "benchmark.csv" "$PDF_NAME";

        [[ -f Rplots.pdf ]] && rm Rplots.pdf
        [[ -f benchmark.csv ]] && rm benchmark.csv
    )
}

accuracy_plot() {
    PDF_NAME="accuracy.pdf";
    (   cd "$SIMULATION_DIR";

        [[ -f $PDF_NAME ]] && rm $PDF_NAME
        cat SK5NHBGXX_measure.json|./simulate.py csv > SK5NHBGXX_measure.csv
        cat SK5NHBGXY_measure.json|./simulate.py csv > SK5NHBGXY_measure.csv

        ./accuracy_plot.R "$PDF_NAME"

        [[ -f Rplots.pdf ]] && rm Rplots.pdf
        [[ -f SK5NHBGXX_measure.csv ]] && rm SK5NHBGXX_measure.csv
        [[ -f SK5NHBGXY_measure.csv ]] && rm SK5NHBGXY_measure.csv
    )
}
