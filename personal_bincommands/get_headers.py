#!/usr/bin/env python

##Script to select columns of interest

###Modules
from argparse import ArgumentParser

#Examples 
#./substract_columns.py  -INPUT coll-pied-merged.gtf  -OUTPUT coll-pied-merged.bed -cols 1,4,5



##Parser and paramters

parser=ArgumentParser(description='To extract the header : all line beginning by a hash in any file')

parser.add_argument("-INPUT",												
                  help="input filename",dest="file_in" )




args = parser.parse_args()
file_input= args.file_in


# say column and  VCF or no


file_in=open(file_input,"r")
file_out=open("headers_"+file_input,"w")




#Read the file
line=file_in.readline()#Read the firsst
while line.startswith("#"):
	file_out.write(line)
	line=file_in.readline()




