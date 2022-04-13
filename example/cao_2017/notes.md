# Analyzing sci-RNA-seq with Pheniqs

This is an analysis of **experiment 3** from the [Cao at el 2017 paper](https://www.science.org/doi/10.1126/science.aam8940) titled *Comprehensive single-cell transcriptional profiling of a multicellular organism*.

## Experiment 3: Single cell RNA profiling of C. elegans

This experiment has two rounds of combinatorial barcoding. The first round adds one of 96 10bp barcodes and an 8bp UMI. The second adds one of 960 20bp Nextera barcodes.

The experiment starts with 150,000 larvae synchronized at the L2 stage dissociated into single-cell suspensions. Cells are split over 6 plates with 96 wells each (a total of 6x96=576 wells). ~1000 C.elegans cells + ~1000 HEK293T cells (as internal controls) are placed in each well. Each of the 96 wells contains a primer with a unique 10bp RT barcode. The set of 96 barcodes is reused in each of the 6 plates. This is the first round of two combinatorial barcoding rounds. The primer also appends an 8bp UMI. In total this operates on 576,000 C.eleganse cells + 576,000 HEK293T cells with one of the 96 RT barcodes and a UMI.

Cells are then pooled and FACS sorted, gating on DNA content to distinguish between C.eleganse and HEK293T cells, and DAPI stain such that singlets are discriminated from doublets. Second strand is synthesized to yield cDNA.

The cells are then sorted into 10 plates with 96 well each for PCR barcoding with dual index Nextera barcodes. 92 out of the 96 wells in each plate contain 140 C.eleganse cells only (a total of 92x10x140 = 128,800 C.eleganse cells). 4 wells in each plate contain a mixture of 140 C.eleganse + 10 HEK293T cells (a total of 4x9x140 = 5,040 C.eleganse cells + 4x9x10 = 360 HEK293T cells), but those 4 wells of plate 4 are not present in the sample sheets.

According to the paper, the experiment yields 42,035 C. elegans single-cell transcriptomes (UMI counts per cell for protein-coding genes ≥ 100). 94% of reads mapped to the expected strand of genic regions (92% exonic and 2% intronic).

## Sample sheet nomenclature

There are a total of 960 i7+i5 combinations. Each is 10bp long resulting in 960 unique 20bp sequences. i7 has a space of 96 combinations and i5 has a space of 48 combinations which yields a total dot product of 4,608 combinations. 960 out of the 4,608 combinations are explicitly selected to uniquely identify each of the 960 wells in the second round of barcoding.

Overall the experiment uses 96 combinations for the first round, RT, barcode and 960 for the second round with a total of 92,160 potential unique cellular barcodes. In the first round cells are split over 6 plates with 96 wells each but use the same 96 RT barcodes which might have to do with efficiency of the reaction but does not actually yield 576 unique combinations.
