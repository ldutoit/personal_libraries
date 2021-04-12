#!/usr/bin/env python2

##Writen by Ludovic Dutoit, April 2021, dutoit.ludovic@gmail.com for any questions

# USAGE EXAMPLE
# python str2vcfPOP.py   INPUT OUTPUT


##LOAD MODULES

import argparse


#parsing argument
parser = argparse.ArgumentParser() # add the parser
parser.add_argument("input",help="input STRUCTURE DATA  file") # add the parser
parser.add_argument("output",help="output VCF file") # add the parser

args = parser.parse_args()




##PROGRAM
print ("WARNING: this script shall never been used to any other purpose than transforming a structure file 2 lines per ind to a vcf it does not encode proper genotype likelihood, qualities and depth information\n\n")
print ("WARNING: Missing data should be encoded as -9\n\n")

## Translation dictionaries for genotypes

dict_individuals = {}

tr={"-9":".","1":"0","2":"1"} # to translate genotypes to vcf alleles
tr_random={"0/1":":19:13,6:40:-14.66,-0.00,-42.03",
			"0/0":":10:10,0:40:-0.00,-4.64,-46.56",
			"1/1":":24:0,24:40:-109.55,-6.60,-0.00",
			"1/0":":19:6,13:40:-14.66,-0.00,-42.03",
			"./.":""
			}


#Load the data ind per ( line by line )
print (" reading "+ args.input+ " ...")
i=0
with open(args.input) as f:
	for line in f:
		i+=1
		if i==1: loci =line.strip().split()
		if i>1:
			info = line.strip().split("\t")
			ind=  info[0].replace(" ","_")
			if not ind in dict_individuals.keys():
				gen = [tr[x] for x in info[1:]]
				dict_individuals[ind]=[gen]
			elif ind in dict_individuals.keys():
				gen = [tr[x] for x in info[1:]]
				dict_individuals[ind].append(gen)

##OUTPUT  the data
print (" writing  "+ args.output+ " ...")

output = open(args.output,"w")


# write the first header
output.write('##fileformat=VCFv4.2\n##fileDate=....\n##source="Random Ludo"\n##INFO=<ID=AD,Number=R,Type=Integer,Description="Total Depth for Each Allele">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##INFO=<ID=loc_strand,Number=1,Type=Character,Description="Genomic strand the corresponding Stacks locus aligns on">\n')

#write the ind headers
output.write("\t".join(["#CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT"])+"\t" + "\t".join(dict_individuals.keys())+"\n")

# write the file locus by locus
for i,loc in enumerate(loci):
	chrom,pos,alleles = loc.split("-")
	newline = "chr%s\t%s\t.\t%s\t%s\t.\tPASS\tNS=75;AF=0.007\tGT" %(chrom,pos,alleles.split("/")[0],alleles.split("/")[1])
	genotypes =[dict_individuals[ind][0][i]+"/"+dict_individuals[ind][1][i] for ind in dict_individuals.keys()]
	linegen="\t".join(genotypes)+"\n"
	newline = newline + "\t"+linegen
	output.write(newline)

output.close()

print ("DONE")
