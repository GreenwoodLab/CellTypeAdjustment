#!/bin/bash

# Need to output "datafileTmp.txt" from the R script "ct_adjustment_example.R" before running CellCDecon
# Instructions on how to install CellCDecon can be found at https://github.com/jameswagner/CellCDecon

./CellCDecon -k 3 -n 46 -f datafileTmp.txt
