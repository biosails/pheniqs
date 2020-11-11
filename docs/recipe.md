---
layout: default
title: "Standard barcode configurations"
permalink: /recipe
id: recipe
---

A library of common barcode sets that are commercially available. You can import those configuration files in your own and extend the decoder definitions. In many cases your expremint will not contain the entire barcode set so you should remove the barcodes that don't apply or alternatively set their `concentration` to `0`.

[IDT for Ilumina Nextera DNA UD]({{ site.github.repository_url }}/blob/master/example/illumina/idt_for_ilumina_nextera_dna_ud.json) Dual indexed with two 10 nucleotide long barcode seqeunces. Segment 0 is the i7 index and segement 1 is the i5 index. Workflow A is used by NovaSeq, MiSeq, HiSeq 2000, HiSeq 2500 i5 is as specified. Workflow B is used by iSeq, MiniSeq, NextSeq, HiSeq 3000,HiSeq 4000 where i5 is reverse complemented.

[IDT for Ilumina TruSeq DNA and RNA UD]({{ site.github.repository_url }}/blob/master/example/illumina/idt_for_ilumina_truseq_dna_and_rna_ud.json) Dual indexed with two 8 nucleotide long barcode seqeunces. Segment 0 is the i7 index and segement 1 is the i5 index. Workflow A is used by NovaSeq, MiSeq, HiSeq 2000, HiSeq 2500 i5 is as specified. Workflow B is used by iSeq, MiniSeq, NextSeq, HiSeq 3000,HiSeq 4000 where i5 is reverse complemented.

[Illumina AmpliSeq]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_ampliseq.json) The i5 (for both Workflow A and Workflow B, with reverse complemented sequence() and i7 are provided separately.

[Illumina Nextera DNA Indexes]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_nextera_dna_indexes.json)

[Illumina TruSeq Amplicon]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_truseq_amplicon.json)

[Illumina TruSeq DNA and RNA CD]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_truseq_dna_and_rna_cd.json)

[Illumina Truseq Small RNA]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_truseq_small_rna.json)

[Illumina Trusight Amplicon]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_trusight_amplicon.json)

[Illumina TruSight DNA Enrichment]({{ site.github.repository_url }}/blob/master/example/illumina/illumina_trusight_dna_enrichment.json)

[TruSeq targeted RNA Expression]({{ site.github.repository_url }}/blob/master/example/illumina/truseq_targeted_rna_expression.json)
