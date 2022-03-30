#!/usr/bin/env python2

## Written by Ludovic Dutoit, improved by Oasis Who
## Note, one line per individual, 2 columns per locus, set POPDATA and LABEL as one in mainparams
#import modules
from __future__ import print_function
import argparse
import sys
import os
#third parties module
import vcf # pip install pyVCF


# this file convert a vcf file into a structure file. require pyVCF library.
## example usage
#./vcf2structure.py --input YOURVCF.vcf --output YOUROUTPUT.str

#FUNCTIONS



def errprint(*args, **kwargs):
    ''' print to stderr not stdout'''
    print(*args, file=sys.stderr, **kwargs)

#parser
parser = argparse.ArgumentParser() # add the parser
parser.add_argument("input",help="input VCF file") # add the parser
parser.add_argument("output",help="output STRUCTURE DATA file") # add the parser

args = parser.parse_args()




dict_alleles = {"0/0":"1\t1","0/1":"1\t2","1/0":"1\t2","1/1":"2\t2","./.":"-9\t-9"}

input_vcf=vcf.Reader(fsock=None, filename=args.input, compressed=False, prepend_chr="False", strict_whitespace=False)#open the vcf parser


list_snps = []
nsites = 0
gen_dict = {ind:[] for ind in input_vcf.samples } 

#store all the genotypes and loci names
for site in input_vcf:
	list_snps.append( site.CHROM+"_"+str(site.POS)) # chr_pos first allele
	for i in range(len(gen_dict.keys())):
		gen_dict[site.samples[i].sample].append(dict_alleles[site.samples[i]["GT"]])

#Write the strcture file
output = open(args.output,"w")
output.write("\t".join(list_snps)+"\n")
for ind in gen_dict.keys():
	#print (gen_dict[ind])
	#print (ind)
	output.write("\t".join([ind]+[str(1)]+gen_dict[ind])+"\n")
output.close()
