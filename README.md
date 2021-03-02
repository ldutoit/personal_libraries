# Personal libraries

This repository contains libraries that I am using regularly in my code. 
Some of it is specific to the UPPMAX cluster system in Uppsala ( might run on other slurm based servers).
For any questions, please email dutoit.ludovic@gmail.com

### stand alone scripts

**swap_min_max.py**  Swap column to have always the first one with the min and the second one with the max.

**vcf2structure.py** convert a vcf to a structure data file

**str2vcfPOP.py** Conversion of a structure file to a vcf populations file for STACKS. Warning: very inadequate for anything else.

**PCR_window_design.py** Find windows suitable for amplification in a fasta alignment. Take into account custom primer design conditions.

**haplo_network.md** Little tutorial to run haplotype networks through R.

### windows_tools.py

Python library that contain tools to work with bedfiles and extract features, include the WindoWBed class to refer to any genomic interval and the Bed files to store windows.


### vcf_tools.py
Python library used to process vcf during my PhD. Some of the main functions are:

**extract_pi_double_vcf_bed()**
Append pi information line by line to a bed file. It uses two vcfs as input ( one for genotypes, and one for depth at all sites)

**count_sites_under_condition_vcf()**
Check the number of sites that respect a condition for a sub part of a vcf file.

**all_freq_spectrum_vcf_bed()**
Append the folded allele frequency spectrum information line by line to a bed file filtering the vcf file according to user based conditions. It uses two vcfs as input ( one for genotypes, and one for depth at all sites)( one for genotypes, and one for depth at all sites)

### afs_stats.py

Calculate Tajima's D, theta Watterson and nucleotide diversity from the allelic frequency spectrum. The main function is **afs_stats**.


### fasta_tools.py

Contains tools to work with fasta format ( parser/writer) as well as assembly information such as scaffold and chromosomes information specific to the flycatcher. This library is dependent on file available in UPPMAX and developped for members of the Ellegren lab.

### R_plotting.R
contains some basic shortcuts I use to plots in R


## Drafts

This folder contains coded skelettons of script that can be useful in the future but that do nothing concrete as is. A good example is a basic argument parser for python that handle stderr and stdout separately.


