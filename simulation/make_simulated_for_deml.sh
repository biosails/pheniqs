#!/usr/bin/env zsh

SAMTOOLS_THREADS=2
SOURCE_DATA="/volume/albireo/waverly/deml/pheniqs/A5KVK.cram"
SIMULATE_PY="/Users/lg/code/biosails/pheniqs/tool/simulate.py"
TODEML_PY="/Users/lg/code/biosails/pheniqs/tool/todeml.py"
SIMULATE_PY_PROFILE="/Users/lg/code/biosails/pheniqs/test/synthetic/simulation.json"
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_24.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_24_for_deml.bam

samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_002408.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_002408_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_010169.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_010169_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_020160.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_020160_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_029793.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_029793_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_039433.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_039433_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_048776.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_048776_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_057784.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_057784_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_066878.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_066878_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_075641.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_075641_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_084253.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_084253_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_092736.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_092736_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_101145.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_101145_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_109136.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_109136_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_117206.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_117206_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_125214.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_125214_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_132844.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_132844_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_140526.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_140526_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_147897.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_147897_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_155148.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_155148_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_162508.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_162508_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_169414.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_169414_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_176525.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_176525_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_183351.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_183351_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_190047.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_190047_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_196708.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_196708_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_203133.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_203133_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_209702.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_209702_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_215989.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_215989_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_222165.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_222165_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_228286.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_228286_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_234293.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_234293_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_240130.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_240130_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_245970.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_245970_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_251751.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_251751_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_257444.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_257444_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_263088.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_263088_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_268502.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_268502_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_273824.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_273824_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_279133.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_279133_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_284395.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_284395_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_289590.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_289590_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_294594.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_294594_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_299673.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_299673_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_304587.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_304587_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_309421.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_309421_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_314206.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_314206_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_318938.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_318938_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_323720.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_323720_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_328343.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_328343_for_deml.bam
samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_332825.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_332825_for_deml.bam
# samtools view --threads $SAMTOOLS_THREADS A5KVK_n50_nu1_337214.bam|$TODEML_PY|samtools view --threads $SAMTOOLS_THREADS -b > A5KVK_n50_nu1_337214_for_deml.bam
