#!/bin/bash

EXTENSION='fa'
KMERSIZE=9
KMERCLUSTER=10
NUMSPACERS=3
EXTRAOPTIONS=
LOGLEVEL=4
while getopts ":E:K:k:l:f:X:" opt; do
    case $opt in
        E)
            EXTENSION=$OPTARG
            ;;
        K)
            KMERSIZE=$OPTARG
            ;;
        k)
            KMERCLUSTER=$OPTARG
            ;;
        l)
            LOGLEVEL=$OPTARG
            ;;
        X)
            EXTRAOPTIONS=$OPTARG
            ;;
        f)
            NUMSPACERS=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

for f in *.${EXTENSION}; do
    echo "processing file $f"
     crass -o crass_out_${f%.${EXTENSION}} -k $KMERCLUSTER -K $KMERSIZE -l $LOGLEVEL -f $NUMSPACERS $EXTRAOPTIONS $f
done
