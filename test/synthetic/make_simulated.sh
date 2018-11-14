#!/usr/bin/env zsh

SOURCE_DATA="/volume/albireo/waverly/deml/pheniqs/A5KVK.cram"
SIMULATE_PY="/Users/lg/code/biosails/pheniqs/tool/simulate.py"
TODEML_PY="/Users/lg/code/biosails/pheniqs/tool/todeml.py"
SIMULATE_PY_PROFILE="/Users/lg/code/biosails/pheniqs/test/synthetic/simulation.json"
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.24 --output A5KVK_n50_nu1_0024.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_0024_model.json

samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.002408 --output A5KVK_n50_nu1_002408.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_002408_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.010169 --output A5KVK_n50_nu1_010169.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_010169_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.020160 --output A5KVK_n50_nu1_020160.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_020160_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.029793 --output A5KVK_n50_nu1_029793.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_029793_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.039433 --output A5KVK_n50_nu1_039433.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_039433_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.048776 --output A5KVK_n50_nu1_048776.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_048776_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.057784 --output A5KVK_n50_nu1_057784.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_057784_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.066878 --output A5KVK_n50_nu1_066878.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_066878_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.075641 --output A5KVK_n50_nu1_075641.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_075641_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.084253 --output A5KVK_n50_nu1_084253.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_084253_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.092736 --output A5KVK_n50_nu1_092736.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_092736_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.101145 --output A5KVK_n50_nu1_101145.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_101145_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.109136 --output A5KVK_n50_nu1_109136.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_109136_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.117206 --output A5KVK_n50_nu1_117206.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_117206_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.125214 --output A5KVK_n50_nu1_125214.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_125214_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.132844 --output A5KVK_n50_nu1_132844.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_132844_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.140526 --output A5KVK_n50_nu1_140526.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_140526_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.147897 --output A5KVK_n50_nu1_147897.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_147897_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.155148 --output A5KVK_n50_nu1_155148.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_155148_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.162508 --output A5KVK_n50_nu1_162508.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_162508_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.169414 --output A5KVK_n50_nu1_169414.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_169414_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.176525 --output A5KVK_n50_nu1_176525.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_176525_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.183351 --output A5KVK_n50_nu1_183351.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_183351_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.190047 --output A5KVK_n50_nu1_190047.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_190047_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.196708 --output A5KVK_n50_nu1_196708.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_196708_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.203133 --output A5KVK_n50_nu1_203133.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_203133_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.209702 --output A5KVK_n50_nu1_209702.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_209702_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.215989 --output A5KVK_n50_nu1_215989.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_215989_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.222165 --output A5KVK_n50_nu1_222165.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_222165_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.228286 --output A5KVK_n50_nu1_228286.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_228286_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.234293 --output A5KVK_n50_nu1_234293.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_234293_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.240130 --output A5KVK_n50_nu1_240130.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_240130_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.245970 --output A5KVK_n50_nu1_245970.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_245970_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.251751 --output A5KVK_n50_nu1_251751.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_251751_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.257444 --output A5KVK_n50_nu1_257444.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_257444_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.263088 --output A5KVK_n50_nu1_263088.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_263088_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.268502 --output A5KVK_n50_nu1_268502.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_268502_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.273824 --output A5KVK_n50_nu1_273824.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_273824_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.279133 --output A5KVK_n50_nu1_279133.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_279133_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.284395 --output A5KVK_n50_nu1_284395.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_284395_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.289590 --output A5KVK_n50_nu1_289590.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_289590_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.294594 --output A5KVK_n50_nu1_294594.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_294594_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.299673 --output A5KVK_n50_nu1_299673.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_299673_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.304587 --output A5KVK_n50_nu1_304587.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_304587_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.309421 --output A5KVK_n50_nu1_309421.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_309421_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.314206 --output A5KVK_n50_nu1_314206.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_314206_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.318938 --output A5KVK_n50_nu1_318938.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_318938_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.323720 --output A5KVK_n50_nu1_323720.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_323720_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.328343 --output A5KVK_n50_nu1_328343.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_328343_model.json
samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.332825 --output A5KVK_n50_nu1_332825.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_332825_model.json
# samtools view $SOURCE_DATA|$SIMULATE_PY -p deml_n50_nu1 -e 0.337214 --output A5KVK_n50_nu1_337214.sam $SIMULATE_PY_PROFILE > A5KVK_n50_nu1_337214_model.json
