#!/bin/bash

#Running FaST-LMM-EWASher.  Install EWASher from http://research.microsoft.com/en-us/downloads/472fe637-7cb9-47d4-a0df-37118760ccd1/
#Need to output "datafileTmp.txt", "phenfileTmp.txt", and "covarfileTmp.txt" from the R script "ct_adjustment_example.R".

#FLE-script.r is provided in the EWASher installation.  Note that you may need to add the following line to the top of FLE-script.r
#for EWASher to run properly: Sys.setenv(FastLmmUseAnyMklLib=1)

#Also, you will need to make an edit in the file "fastlmm-ewasher.r".  On line 230/231 (in function "CreateFastlmmInputFiles") replace:
# write.table(scores[, 1:100], file = 'pc.txt', quote = FALSE, sep = '\t',
#        row.names = FALSE, col.names = FALSE)
#with:
# write.table(scores, file = 'pc.txt', quote = FALSE, sep = '\t',
#        row.names = FALSE, col.names = FALSE)

Rscript FLE-script.r datafileTmp.txt phenfileTmp.txt -covar covarfileTmp.txt

