#!/usr/bin/env zsh

# 1. dataset in SAM/BAM/CRAM format
# 2. simulate_barcode: use the priors and barcode words to simulate barcode indices in the data
#    with perfect barcodes and the existing quality scores. determine the global error rate on the
#    cycles containing the barcode indice. correct barcode index is encoded in the read id.
# 3. simulate_error: using the global error rate and a desired error rate, simulate errors according
#    to the quality scores after recalibrating them for the new global rate.
# 4. todeml: create a transcoded version of the data in a format deML excepts
# 5. analyze: collect FP/FN/TP/TN statistics and compile a report
# 6. collect: update the global database with results from the experiment
# 7. summarize: extract data in csv according to some transformation

PHENIQS_HOME="$(cd ~/code/biosails/pheniqs && pwd)"
ORIGINAL_DATA="$(cd ~/Downloads && pwd)/A5KVK_tiny.bam"
echo "PHENIQS_HOME: $PHENIQS_HOME"
echo "ORIGINAL_DATA: $ORIGINAL_DATA"

simulate_barcode() {
# $1 : input SAM
# $2 : profile
samtools view $1 | \
$PHENIQS_HOME/tool/benchmark.py simulate_barcode \
-c $PHENIQS_HOME/simulation/core.json -p $2 \
-r $2_simulated_barcode_model.json | \
samtools view -b > $2_simulated_barcode.bam
}

simulate_error() {
# $1 : profile
# $2 : desired error rate
samtools view $1_simulated_barcode.bam | \
$PHENIQS_HOME/tool/benchmark.py simulate_error -e 0.$2 \
-m $1_simulated_barcode_model.json \
-r $1_$2_simulated_error_model.json | \
samtools view -b > $1_$2_simulated_error.bam
}

todeml() {
# $1 : profile
# $2 : desired error rate
samtools view $1_$2_simulated_error.bam | \
$PHENIQS_HOME/tool/benchmark.py todeml | \
samtools view -b > $1_$2_simulated_error_deml.bam
}

sense_prior() {
# $1 : profile
# $2 : desired error rate
$PHENIQS_HOME/tool/benchmark.py sense_prior \
-c pheniqs_annotate.json \
-i $1_$2_simulated_error.bam \
-r $1_$2_prior_estimate.json
}

adjust_prior() {
# $1 : profile
# $2 : desired error rate
$PHENIQS_HOME/tool/benchmark.py adjust_prior \
-c pheniqs_annotate.json \
-r $1_$2_prior_estimate.json \
-m $1_$2_simulated_error_model.json \
-a $1_$2_adjusted_annotate.json
}

demux_pheniqs() {
# $1 : profile
# $2 : desired error rate
$PHENIQS_HOME/tool/benchmark.py demux_pheniqs \
-c $1_$2_adjusted_annotate.json \
-i $1_$2_simulated_error.bam \
-o $1_$2_simulated_pheniqs_demux.bam \
-r $1_$2_pheniqs_demux_report.json
}

demux_deml() {
# $1 : profile
# $2 : desired error rate
$PHENIQS_HOME/tool/benchmark.py demux_deml \
-c deml_index.txt \
-i $1_$2_simulated_error_deml.bam \
-o $1_$2_simulated_deml_demux.bam \
-r $1_$2_deml_summary.txt
}

analyze_pheniqs() {
# $1 : profile
# $2 : desired error rate
samtools view $1_$2_simulated_pheniqs_demux.bam | \
$PHENIQS_HOME/tool/benchmark.py analyze_pheniqs \
-m $1_$2_simulated_error_model.json \
-r $1_$2_simulated_pheniqs_analysis.json
}

analyze_deml() {
# $1 : profile
# $2 : desired error rate
samtools view $1_$2_simulated_deml_demux.bam | \
$PHENIQS_HOME/tool/benchmark.py analyze_deml \
-m $1_$2_simulated_error_model.json \
-r $1_$2_simulated_deml_analysis.json
}

execute_experiment() {
  simulate_error $1 $2
  todeml $1 $2
  sense_prior $1 $2
  adjust_prior $1 $2
  demux_pheniqs $1 $2
  demux_deml $1 $2
  analyze_pheniqs $1 $2
  analyze_deml $1 $2
}

simulate_barcode $ORIGINAL_DATA A5KVK

execute_experiment A5KVK 0001
execute_experiment A5KVK 002408
execute_experiment A5KVK 010169
execute_experiment A5KVK 020160
execute_experiment A5KVK 029793
execute_experiment A5KVK 039433
execute_experiment A5KVK 048776
execute_experiment A5KVK 057784
execute_experiment A5KVK 066878
execute_experiment A5KVK 075641
execute_experiment A5KVK 084253
execute_experiment A5KVK 092736
execute_experiment A5KVK 101145
execute_experiment A5KVK 109136
execute_experiment A5KVK 117206
execute_experiment A5KVK 125214
execute_experiment A5KVK 132844
execute_experiment A5KVK 140526
execute_experiment A5KVK 147897
execute_experiment A5KVK 155148
execute_experiment A5KVK 162508
execute_experiment A5KVK 169414
execute_experiment A5KVK 176525
execute_experiment A5KVK 183351
execute_experiment A5KVK 190047
execute_experiment A5KVK 196708
execute_experiment A5KVK 203133
execute_experiment A5KVK 209702
execute_experiment A5KVK 215989
execute_experiment A5KVK 222165
execute_experiment A5KVK 228286
execute_experiment A5KVK 234293
execute_experiment A5KVK 240130
execute_experiment A5KVK 245970
execute_experiment A5KVK 251751
execute_experiment A5KVK 257444
execute_experiment A5KVK 263088
execute_experiment A5KVK 268502
execute_experiment A5KVK 273824
execute_experiment A5KVK 279133
execute_experiment A5KVK 284395
execute_experiment A5KVK 289590
execute_experiment A5KVK 294594
execute_experiment A5KVK 299673
execute_experiment A5KVK 304587
execute_experiment A5KVK 309421
execute_experiment A5KVK 314206
execute_experiment A5KVK 318938
execute_experiment A5KVK 323720
execute_experiment A5KVK 328343
execute_experiment A5KVK 332825
execute_experiment A5KVK 337214
