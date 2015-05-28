#!/usr/bin/env python

##Script to select columns of interest

###Modules
from argparse import ArgumentParser

#Examples 
#./swap_min_max.py  -INPUT fAlb15.RM.chrompos.bed  -OUTPUT fAlb15.RM.chrompos.sorted.bed -cols 2,3



##Parser and paramters

parser=ArgumentParser(description='To select the columns of interest only in any kind of tab separated  files. If start by double-dashed lines, they are reported in the output without columns substractions(i2. allow to parse vcf.')

parser.add_argument("-INPUT",												
                  help="input filename",dest="file_in" )
parser.add_argument("-OUTPUT",												
                  help="Ooutput filename",dest="file_out")
parser.add_argument("-cols",												
                  help="a comma separated list of the columns to swap",dest="col_number" ,type=str)
 




args = parser.parse_args()
file_input= args.file_in
file_output = args.file_out
col_number = args.col_number


# say column and  VCF or not


#Parameters
#file_input="coll-pied-merged.gtf"#input filename
#file_output="coll-pied-merged.bed"#output filename
#col_number="2,3"#A string, comma separated with col number (i.e "1,2,3,5,"")column 1 is called with 1 and not with zero





file_in=open(file_input,"r")
file_out=open(file_output,"w")


###Transform the string of columns into  python position list
col_number_not_shifted=col_number.split(",")
col_number_list=[]
for i in col_number_not_shifted:
	col_number_list.append(str(int(i)-1))


#Read the file


line=file_in.readline()#Read the firsst
while line.startswith("##"):
	oldline=line#just to retrieve
	file_output.write(oldline)
	line=file_in.readline()

#Read the rest of the file
while len(line)>0:
	to_swap=line.split()
	swap=to_swap
	minimum=str(min(int(to_swap[int(col_number_list[0])]),int(to_swap[int(col_number_list[1])])))
	maximum=str(max(int(to_swap[int(col_number_list[0])]),int(to_swap[int(col_number_list[1])]))	)
	swap[int(col_number_list[0])] =minimum
	swap[int(col_number_list[1])]=maximum
	to_write="\t".join(swap)+"\n"
	file_out.write(to_write)
	line=file_in.readline()

file_out.close()



