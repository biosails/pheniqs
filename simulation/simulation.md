Decoder:
    Noise:
        TP : ['barcode'][0]['noise'][qc]['TP']
             a noise read was correctly classified as noise

        FN : ['barcode'][0]['noise'][qc]['FN']
             a noise read was incorrectly classified to a real barcode

        FP : ['barcode'][0]['noise'][qc]['FP']
             a real read was classified as noise

Real barcode accumulator:
  TP : real read is correctly classified
  FN :
      real: a read from this barcode was incorrectly classified to another real barcode
      noise: a read from this barcode was incorrectly classified as noise
  FP :
      real: a real read from another barcode was incorrectly classified to this barcode
      noise: a noise read was incorrectly classified to this barcode
