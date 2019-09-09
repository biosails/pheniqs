#!/bin/zsh

PROJECT_HOME="$HOME/code/private/pheniqs";
SIMULATION_DIR="$PROJECT_HOME/benchmark/HK5NHBGXX/simulation"
BENCHMAKR_DIR="$PROJECT_HOME/benchmark";
MANUSCRIPT_DIR="$PROJECT_HOME/manuscript/bioinformatics";
TEMPLATE_DIR="$PROJECT_HOME/manuscript/template/bioinformatics"

benchmark_plot() {
    PDF_NAME="benchmark.pdf";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$BENCHMAKR_DIR";
        cat benchmark.log | ./benchmark.py > benchmark.csv;
        ./benchmark.R "benchmark.csv" "$PDF_NAME";

        cp $PDF_NAME $MANUSCRIPT_DIR/$PDF_NAME
        cp $PDF_NAME $TEMPLATE_DIR/$PDF_NAME

        [[ -f $PDF_NAME ]] && rm $PDF_NAME
        [[ -f Rplots.pdf ]] && rm Rplots.pdf
        [[ -f benchmark.csv ]] && rm benchmark.csv
    )
}

accuracy_plot() {
    PDF_NAME="accuracy.pdf";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$SIMULATION_DIR";

        [[ -f $PDF_NAME ]] && rm $PDF_NAME

        cat SK5NHBGXX_measure.json|./simulate.py csv > SK5NHBGXX_measure.csv
        cat SK5NHBGXY_measure.json|./simulate.py csv > SK5NHBGXY_measure.csv
        ./accuracy_plot.R "$PDF_NAME"

        cp $PDF_NAME $MANUSCRIPT_DIR/$PDF_NAME
        cp $PDF_NAME $TEMPLATE_DIR/$PDF_NAME

        # clean up
        [[ -f $PDF_NAME ]] && rm $PDF_NAME
        [[ -f Rplots.pdf ]] && rm Rplots.pdf
        [[ -f SK5NHBGXX_measure.csv ]] && rm SK5NHBGXX_measure.csv
        [[ -f SK5NHBGXY_measure.csv ]] && rm SK5NHBGXY_measure.csv
    )
}

biological_plot() {
    BIOLOGICAL_DATA="/volume/albireo/waverly/HK5NHBGXX/HK5NHBGXX_R.csv.bz2";
    PDF_NAME="biological.pdf";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$BENCHMAKR_DIR/HK5NHBGXX";

        [[ -f $PDF_NAME ]] && rm $PDF_NAME

        cat error.json | ./error.py > error.tsv;
        ./analyze_sam.R $BIOLOGICAL_DATA error.tsv $PDF_NAME

        cp $PDF_NAME $MANUSCRIPT_DIR/$PDF_NAME
        cp $PDF_NAME $TEMPLATE_DIR/$PDF_NAME

        # clean up
        [[ -f $PDF_NAME ]] && rm $PDF_NAME
        [[ -f Rplots.pdf ]] && rm Rplots.pdf
        [[ -f error.csv ]] && rm error.csv
    )
}

manuscript_bioinformatics() {
    PDF_NAME="pheniqs.pdf";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$TEMPLATE_DIR";
        [[ -f pheniqs.aux ]] && rm pheniqs.aux;
        [[ -f pheniqs.log ]] && rm pheniqs.log;
        [[ -f pheniqs.bbl ]] && rm pheniqs.bbl;
        [[ -f pheniqs.blg ]] && rm pheniqs.blg;
        [[ -f pheniqs.pdf ]] && rm pheniqs.pdf;

        pdflatex pheniqs && bibtex8 pheniqs && pdflatex pheniqs && pdflatex pheniqs;
        [[ -f pheniqs.pdf ]] && mv pheniqs.pdf "$MANUSCRIPT_DIR/$PDF_NAME"

        [[ -f pheniqs.aux ]] && rm pheniqs.aux;
        [[ -f pheniqs.log ]] && rm pheniqs.log;
        # [[ -f pheniqs.bbl ]] && rm pheniqs.bbl;
        [[ -f pheniqs.blg ]] && rm pheniqs.blg;
    )
}

supplementary_bioinformatics() {
    PDF_NAME="supplementary.pdf"
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$TEMPLATE_DIR";
        [[ -f supplementary.aux ]] && rm supplementary.aux;
        [[ -f supplementary.log ]] && rm supplementary.log;
        [[ -f supplementary.pdf ]] && rm supplementary.pdf;

        pdflatex supplementary;
        [[ -f supplementary.pdf ]] && mv supplementary.pdf "$MANUSCRIPT_DIR/$PDF_NAME"

        [[ -f supplementary.aux ]] && rm supplementary.aux;
        [[ -f supplementary.log ]] && rm supplementary.log;
    )
}

latex_bioinformatics() {
    LATEX_NAME="pheniqs.tex";
    BBL_NAME="pheniqs.bbl";
    REF_NAME="reference.bib";
    BIOINFO_NAME="bioinfo.cls";
    [[ -f "$MANUSCRIPT_DIR/$LATEX_NAME" ]] && rm "$MANUSCRIPT_DIR/$LATEX_NAME"
    [[ -f "$MANUSCRIPT_DIR/$BBL_NAME" ]] && rm "$MANUSCRIPT_DIR/$BBL_NAME"
    [[ -f "$MANUSCRIPT_DIR/$REF_NAME" ]] && rm "$MANUSCRIPT_DIR/$REF_NAME"
    [[ -f "$MANUSCRIPT_DIR/$BIOINFO_NAME" ]] && rm "$MANUSCRIPT_DIR/$BIOINFO_NAME"
    (   cd "$TEMPLATE_DIR";
        [[ -f $LATEX_NAME ]] && cp $LATEX_NAME "$MANUSCRIPT_DIR/$LATEX_NAME";
        [[ -f $REF_NAME ]] && cp $REF_NAME "$MANUSCRIPT_DIR/$REF_NAME";
        [[ -f $BBL_NAME ]] && cp $BBL_NAME "$MANUSCRIPT_DIR/$BBL_NAME";
        [[ -f $BIOINFO_NAME ]] && cp $BIOINFO_NAME "$MANUSCRIPT_DIR/$BIOINFO_NAME";
    )
}

manuscript_biorxiv() {
    PDF_NAME="pheniqs.pdf";
    TEMPLATE_DIR="$PROJECT_HOME/manuscript/template/biorxiv"
    MANUSCRIPT_DIR="$PROJECT_HOME/manuscript/biorxiv";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$TEMPLATE_DIR";
        [[ -f pheniqs.aux ]] && rm pheniqs.aux;
        [[ -f pheniqs.log ]] && rm pheniqs.log;
        [[ -f pheniqs.bbl ]] && rm pheniqs.bbl;
        [[ -f pheniqs.blg ]] && rm pheniqs.blg;
        [[ -f pheniqs.pdf ]] && rm pheniqs.pdf;

        pdflatex pheniqs && biber pheniqs && pdflatex pheniqs && pdflatex pheniqs;
        [[ -f pheniqs.pdf ]] && mv pheniqs.pdf "$MANUSCRIPT_DIR/$PDF_NAME"

        [[ -f pheniqs.aux ]] && rm pheniqs.aux;
        [[ -f pheniqs.log ]] && rm pheniqs.log;
        [[ -f pheniqs.bbl ]] && rm pheniqs.bbl;
        [[ -f pheniqs.blg ]] && rm pheniqs.blg;
    )
}

supplementary_biorxiv() {
    PDF_NAME="supplementary.pdf"
    TEMPLATE_DIR="$PROJECT_HOME/manuscript/template/biorxiv"
    MANUSCRIPT_DIR="$PROJECT_HOME/manuscript/biorxiv";
    [[ -f "$MANUSCRIPT_DIR/$PDF_NAME" ]] && rm "$MANUSCRIPT_DIR/$PDF_NAME"
    (   cd "$TEMPLATE_DIR";
        [[ -f supplementary.aux ]] && rm supplementary.aux;
        [[ -f supplementary.log ]] && rm supplementary.log;
        [[ -f supplementary.pdf ]] && rm supplementary.pdf;

        pdflatex supplementary;
        [[ -f supplementary.pdf ]] && mv supplementary.pdf "$MANUSCRIPT_DIR/$PDF_NAME"

        [[ -f supplementary.aux ]] && rm supplementary.aux;
        [[ -f supplementary.log ]] && rm supplementary.log;
    )
}

clean_pdflatex() {
    [[ -f pheniqs.bcf ]] && rm pheniqs.bcf
    [[ -f pheniqs.aux ]] && rm pheniqs.aux
    [[ -f pheniqs.log ]] && rm pheniqs.log
    [[ -f pheniqs.run.xml ]] && rm pheniqs.run.xml
    [[ -f pheniqs.bbl ]] && rm pheniqs.bbl
    [[ -f pheniqs.blg ]] && rm pheniqs.blg
    [[ -f pheniqs.bcfdoc ]] && rm pheniqs.bcf
}

make_pdflatex() {
    [[ -f pheniqs.pdf ]] && rm pheniqs.pdf

    pdflatex pheniqs
    pdflatex pheniqs
    bibtex pheniqs
    pdflatex pheniqs
    pdflatex pheniqs
}
