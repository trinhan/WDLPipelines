#!/bin/bash

## Version 1.1m 28-1-2023
## This is a bash script to perform GCP operations:
## OPTION move:
##         bash file to transfer gcp files from location 1 to location 2
##         required: a tab-delim table with original location - new location
## OPTION delete:
##        bash file to delete samples
##         How to use:
## Pre-requisites:
## 1. Must have gsutil installed and correct permissions for the bucekt of interest

## Enter usage and the required variables
usage()
{
    echo "usage: 
           bash bash_gcp_reorganise.sh 
            -i file list of interest [ either a 2 column table from moveFiles.py or 1 column list from ListUntabled.py  ]
            -m run mode [can either be 'move' or 'delete']
            -t test mode [only prints the outputs to screen, does not actually perform the operations]
"
}

while getopts i:m:t flag
do
    case "${flag}" in
        i) i=${OPTARG};;
        m) m=${OPTARG};;
        t) s=0; echo "Running test mode, just print to screen";;
    esac
done

## Set up exit errors if the above files do not exist
[[ -z "$m" ]] && { usage; echo "No run mode selected"; exit 1; }
[[ -z "$i" ]] && { usage; echo "No input file give"; exit 1; }
[[ -z "$s" ]] && { echo "Performing actual operations"; s=1; }


## Set up green text colour to highlight inputs
GREEN='\033[0;32m'
NC='\033[0m'

## To do: check the run mode is consistent with the set up
## ie. a 2 column table for move and 1 column for delete
## Also check that gsutil is installed

## Display the parameters as a check before proceeding correct:
echo -e "Proceed to ${GREEN}$m${NC} all the files in attached file ${GREEN}$i${NC}?" 
read -p "Note this cannot be undone!" yn
case $yn in
    [Yy]* ) echo 'Proceeding with operations'
        while read -r -a line; do
            case $m in 
                'move' )
                    if [ -z ${line[1]} ]
                    then
                        :
                    else
                        echo "moving sample ${line[0]} to ${line[1]}"
                        if [ $s -eq "1" ]
                        then
                            #echo 'run'
                            gsutil mv "${line[0]}" "${line[1]}"
                        fi
                    fi;;
                'delete' )
                    echo "deleting sample ${line[0]}"
                    if [ $s -eq "1" ]
                    then
                        echo run
                        gsutil rm "${line[0]}"
                    fi;;
            esac
        done < $i ;;
    [Nn]* ) echo 'Halting operations'; exit 1;;
        * ) echo "Please answer y or n.";;
esac


