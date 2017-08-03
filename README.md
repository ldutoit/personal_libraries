# Personal libraries

This repository contain library that I am using regularly in my code. 
Some of it is specific to the UPPMAX cluster system in Uppsala.
For any questions, please email dutoit.ludovic@gmail.com


## windows_tools.py

Python library that contain tools to work with bedfiles and extract features, include the WindoWBed class to refer to any genomic interval and the Bed files to store windows.


## vcf_tools.py
Python library used to process vcf during my PhD. Some of the main functions:

*extract_pi_double_vcf_bed*
Append pi information line by line to a bed file. It uses two vcfs as input ( one for genotypes, and one for depth at all sites)

*count_sites_under_condition_vcf*
Check the number of sites that respect a condition for a sub part of a vcf file.

*all_freq_spectrum_vcf_bed*
Append the folded allele frequency spectrum information line by line to a bed file filtering the vcf file according to user based conditions. It uses two vcfs as input ( one for genotypes, and one for depth at all sites)( one for genotypes, and one for depth at all sites)



## afs_stats.py

Calculate Tajima's D, theta Watterson and nucleotide diversity from the allelic frequency spectrum. The main function is *afs_stats.py*.

## personal_bincommands
Contains small bin commands. help is always available woth the "-h flag"

sortbig_bed_vcf.sh ### sort a big bed or vcf file
faSplit ## split a fasta file into small files ( developped by UCSC)
maf_convert.py # convert multple alignemtne format maf and psl ( developped by UCSC)
swap_min_max.py # used to swap column to have always the first one with the min and the second one with the max


## fasta_tools.py

contains tools to work with fasta_tools as well as assembly information such as scaffold_linkages for the flycatcher. This library is dependent on file available in UPPMAX and developped for members of the ellegren lab.

## R_plotting.R
contains some basic shortcuts I use to plots in R

## folder_tools.py

Contains a small function to  check logs files out of UPPMAX



