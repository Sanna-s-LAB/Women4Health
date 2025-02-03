#!/bin/bash

# This script run the independence tests between PCoA and categorical variables
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 03/02/2025

for correction in T F #run the script "independence_tests.R" for correction = T and F
do
	echo "Doing with correction = ${correction}"
	Rscript "independence_tests.R" $correction
done

