#!/bin/sh
file=$1
##estimate the insert size from the sam file given as first argument
###Example usage:
estimate_insert_size myfile.sam
head -10000 $1 | awk '{if ($9 > 0) {S+=$9; T+=1}}END{print "Mean: " S/T}'


