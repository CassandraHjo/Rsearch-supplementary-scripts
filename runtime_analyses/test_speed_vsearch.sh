#!/bin/bash

echo "elapsed_seconds,max_rss_kb" > results.txt

FORMAT="%e,%M"

for i in {1..5}; do
    gtime -f "$FORMAT" -a -o results.txt bash vsearch_unoise.sh # Change the bash script for each test
done

