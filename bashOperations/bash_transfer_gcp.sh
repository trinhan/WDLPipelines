#!/bin/bash

## bash file to transfer gcp files from location 1 to location 2
## required: a tab-delim table with original location - new location
## How to use:
## bash bash_transfer_gcp.sh move_samples.txt

while read -r -a line; do
	echo 'sample 1'
	echo "${line[0]}"
	echo 'sample 2'
	echo "${line[1]}"
	gsutil mv "${line[0]}" "${line[1]}"
done < $1