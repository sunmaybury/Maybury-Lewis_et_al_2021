#!/bin/bash

# ================================================================================================
# This script takes the output from GREAT analysis and cleans it up with only gene symbols
# Start with clean text files from running gene.assignment.proximal.distal.1.R in same directory
# ================================================================================================

cd ~/Dropbox/Desktop/ATAC-seq-PEAKS/AvsQ_analysis_2019_v2/GREAT_assign


# =================
# AvsQ.DBsites
# =================

cut -f7 "AvsQ.DBsites.proximal+distal.clean.txt" > "AvsQ.DBsites.proximal+distal.tmp.txt"
sed 's/(.*//' "AvsQ.DBsites.proximal+distal.tmp.txt" | tr -d '[[:blank:]]' > "AvsQ.DBsites.proximal+distal.genes.txt"

cut -f7 "AvsQ.DBsites.proximal.clean.txt" > "AvsQ.DBsites.proximal.tmp.txt"
sed 's/(.*//' "AvsQ.DBsites.proximal.tmp.txt" | tr -d '[[:blank:]]' > "AvsQ.DBsites.proximal.genes.txt"

# =================
# AvsQ.truStable
# =================

cut -f5 "AvsQ.truStable.proximal+distal.clean.txt" > "AvsQ.truStable.proximal+distal.tmp.txt"
sed 's/(.*//' "AvsQ.truStable.proximal+distal.tmp.txt" | tr -d '[[:blank:]]' > "AvsQ.truStable.proximal+distal.genes.txt"

cut -f5 "AvsQ.truStable.proximal.clean.txt" > "AvsQ.truStable.proximal.tmp.txt"
sed 's/(.*//' "AvsQ.truStable.proximal.tmp.txt" | tr -d '[[:blank:]]' > "AvsQ.truStable.proximal.genes.txt"


rm *.tmp.txt

