#!/bin/bash

echo "elapsed_seconds,max_rss_kb" > results.txt

FORMAT="%e,%M"

for i in {1..5}; do
    gtime -f "$FORMAT" -a -o results.txt Rscript rsearch_unoise.R # Change the R script for each test
done

