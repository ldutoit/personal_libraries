#!/usr/bin/env python2
##WRITTEN BY LUDOVIC DUTOIT, AUGUST 2018
#This file aims at detecting windows that are suitable for amplification in a Fasta alignment
#It is full of conditions and the proposed primers should always checked visually in the alignbment using a simple alignment visualisation software


####EXAMPLE USAGE
#python PCR_window_design.py -length_window 310 -length_primer 20  INPUT.fasta summary.txt  OUTPUT.fasta 6 10
# this will look for 310 bp region with 20bp primers in the file INPUT.fasta creating the file OUTPUT.fasta.
#Primers are required to have at least 6 invariable sites at their 3' end and a maximum of 10 degenerated sites each

#usage: PCR_window_design.py [-h]
#                            input_file output_summary output_alignment
#                            length_window length_primer n_inv_primer
#                            n_deg_sites#

#positional arguments:
#  input_file        MULTIFASTA input file
#  output_summary    File with results
#  output_alignment  File with new alignmment( just trimming the end where
#                    there is no info)
#  length_window     length of desired sequence to amplify
#  length_primer     length of primers
#  n_inv_primer      number of sites that are required invariable at the 3' end
#                    of both primers
#  n_deg_sites       number of sites that are allowed to be variable in the
#                    primer#

#optional arguments:
#  -h, --help        show this help message and exit#


#parser
import argparse
parser = argparse.ArgumentParser() # add the parser
parser.add_argument("input_file",help="MULTIFASTA input file") # add the parser
parser.add_argument("output_summary",help="File with results") # add the parser
parser.add_argument("output_alignment",help="File with new alignmment( just trimming the end where there is no info)") # add the parser
parser.add_argument("n_inv_primer",help="number of sites that are required invariable at the 3' end of both primers",type=int) # add the parser
parser.add_argument("n_deg_sites",help="number of sites that are allowed to be variable in the primer",type=int) # add the parser
parser.add_argument("-length_window",help="length of desired sequence to amplify (default: 300)",type=int,default=300) # add the parser
parser.add_argument("-length_primer",help="length of primers (default: 20)",type=int,default=20) # add the parser

args = parser.parse_args()


# classes
class window(object):  # come from_vcf_tools in my personal_libraries
	""" a class to define genomic interval, each object has to be defined by at least a seq and two of the three following parameters: start,end,len. 
	The third is inferred
	If format_input!=BED, it assumes a 1 based inclusive format and transfor it into bed
	>>> WindowBed(seq="chr1",start=0,end=10)
	'chr1    0    10'
	>>> WindowBed(seq="chr1",start=0,length=10)
	'chr1    0    10'
	>>> WindowBed(seq="chr1",start=1,end=10,format_input="non_BED")
	'chr1    0    10'
	>>> WindowBed(seq="chr1",start=1,length=10,format_input="non_BED")
	'chr1    0    10'
	""" 
	def __init__ (self,seq,start,end=None,length=None,format_input="BED",strand=None):
		self.seq=seq
		self.start=int(start)
		self.end=end
		self.length=length
		self.strand=strand
		if  all(v is None for v in [self.end,self.length]):
			raise Exception("Window definition error: end and length have to be defined, provide at least one of them")
		if self.end==None:
			self.end=int(self.start)+int(self.length)
		elif self.length==None:
				self.length=int(self.end)-int(self.start)
		if self.end<self.start: raise Exception("Error: Window end position is before start "+self.seq+" "+str(self.start)+" "+str(self.end))
		if self.length==0:  Warning("Error: Window has length 0 "+self.seq+" "+str(self.start)+" "+str(self.end))
		if  not strand in ["+","-",None]: raise Exception("Error: unauthorized strand: %s" % (str(self.strand)))
	def __repr__(self):
		return repr("    ".join([str(self.seq),str(self.start),str(self.end)]))

IUPAC_alphabet = 	{
  "R" : set(["G","A"]),
  "Y" : set(["T","C"]),
  "S" : set(["C","G"]),
  "W" : set(["T","A"]),
  "K" : set(["T","G"]),
  "M" : set(["C","A"]),
  "B"	:set(["C", "G", "T"]),
  "D"	:set(["A", "G", "T"]),	
  "H"	:set(["A", "C", "T"]),	
  "V"	:set(["A", "C", "G"])}




#Functions

#remove starting "-"
def getpos(pos,dic_multifasta):
	''' return a set with all the base found at position pos in the values of dic_multifasta'''
	bases = set([x[pos] for x in dict_fasta.values()])
	while any(base not in "ATGCN-" for base  in bases): #DEAL WITH IUPAC
		for base in bases:
			if base not in "ATGCN-":
				bases = bases.intersection(IUPAC_alphabet[base])
	assert all([base in "ATGCN-"  for base in bases ])
	return bases

def number_of_variable_site(start,end,dic_multifasta):
	''' return a list with the number of bases found for each position in the window defined from start to end (0 base included), ignores N or - '''
	numb_variable_per_site=[]
	for i in range(start,end+1):
		numb_variable_per_site.append(len([x for x in getpos(i,dic_multifasta) if x in "ATGC"]))
		#print [x for x in getpos(i,dic_multifasta) if x !="-" and x!="N"]
	#print("checked",len(numb_variable_per_site),"sites for variability")
	return numb_variable_per_site

def count_variable_sites(start,end,dic_multifasta):
	''' return the number of variable sites in seq from start to end (0 base included), ignores N or -'''
	n_variable_sites =  len([x for x in number_of_variable_site(start,end,dic_multifasta) if x >1])
	#print(n_variable_sites,"n_variable_sites")
	return len([x for x in number_of_variable_site(start,end,dic_multifasta) if x >1])


def getseqprimer(start,end,dic_multifasta):
	for value in dict_fasta.values(): 
		seq = list(value[start:end])
		if not "-" in value:
			break	
	assert len(seq)==args.length_primer
	num_variable = number_of_variable_site(start,end-1,dic_multifasta)
	for index,item in enumerate(num_variable):
		if item>1: seq[index]="N"
	while "".join(seq).find("-")>-1: # wjhile it is there:
		 seq["".join(seq).find("-")]= "N"
	return "".join(seq)

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',"N":"N"} 
def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

#Read the multi fasta
f = open(args.input_file,'r')
dict_fasta={}
name = ''
for line in f:
    #if your line starts with a > then it is the name of the following sequence
    if line.startswith('>'):
        name = line[1:-1]
        dict_fasta[name]=''
        continue #this means skips to the next line
    #This code is only executed if it is a sequence of bases and not a name.
    dict_fasta[name]+=line.strip()

#check that all the sequences in the fasta are the same length
assert len(set([len(x) for x in dict_fasta.values()]))==1 #
len_sequence = len(dict_fasta.values()[0])
print("The total length of the sequence is "+str(len_sequence)+" bp")


###trim the beginning

i=-1
bases=set(["-"])
while bases == set(["-"]):
	i+=1
	bases = getpos(i,dict_fasta)
	#print(i,bases)

dict_fasta={key:dict_fasta[key][i:]for key in dict_fasta.keys()}
len_sequence = len(dict_fasta.values()[0])

###Trim the end
i=len_sequence
bases=set(["-"])
while bases == set(["-"]):
	i-=1
	bases = getpos(i,dict_fasta)
	#print(i,bases)

dict_fasta={key:dict_fasta[key][:1440]for key in dict_fasta.keys()}
len_sequence = len(dict_fasta.values()[0])
print("The total length after timming of start and end '-' sequence is "+str(len_sequence)+" bp")

#for window by window from the length of the primer, go through the file
windows=[]
for i in range(args.length_primer,len_sequence-args.length_primer-args.length_window):
	print ("start at", i)
	currentwindow = window(seq="",start=i,end=i+args.length_window) 
	currentwindow.primeruppos=[i-args.length_primer,i]
	currentwindow.primerdownpos=[i+args.length_window,i+args.length_window+args.length_primer]
	currentwindow.CHECK_n_inv_primer_up = all(x==1 for x in number_of_variable_site(currentwindow.primeruppos[1]-args.n_inv_primer+1,currentwindow.primeruppos[1],dict_fasta))
	currentwindow.CHECK_n_inv_primer_down =  all(x==1 for x in number_of_variable_site(currentwindow.primerdownpos[0],currentwindow.primerdownpos[0]+args.n_inv_primer-1,dict_fasta))
	currentwindow.n_variable_sites_in_window = count_variable_sites(currentwindow.start,currentwindow.end-1,dict_fasta)
	currentwindow.primerup =getseqprimer(currentwindow.primeruppos[0],currentwindow.primeruppos[1],dict_fasta)
	currentwindow.primerdown =reverse_complement( getseqprimer(currentwindow.primerdownpos[0],currentwindow.primerdownpos[1],dict_fasta)	)
	currentwindow.ndeg_sites_up = currentwindow.primerup.count("N")
	currentwindow.ndeg_sites_down = currentwindow.primerdown.count("N")
	currentwindow.validation = currentwindow.CHECK_n_inv_primer_up and currentwindow.CHECK_n_inv_primer_down and currentwindow.ndeg_sites_down <= args.n_deg_sites and currentwindow.ndeg_sites_up  <= args.n_deg_sites
	windows.append(currentwindow)


#Write window by window ( OK first)
output=open(args.output_summary,"w")
output.write("#We selected windows of %i bp  \n#Primers of %i bp with at least %i invariable sites at the end and max %i degenerated sites in each primer\n#windows passing criteria ( if any) at the beginning of toqhe file\n" % (args.length_window,args.length_primer,args.n_inv_primer,args.n_deg_sites))
output.write("\t".join(["respect_criteria","start","end","n_inv_sites_end_primersup","n_inv_sites_end_primers_down","total_ndeg_primersup","total_ndeg_primersudown","total_variable_sites_in_sequence","primerup","primerdown"])+"\n")
n_ok_windows=0
with open(args.output_summary) as f:
	for window in windows:
		if window.validation:
			output.write("\t".join([str(x) for x in [window.validation,window.start+1,window.end+1,window.CHECK_n_inv_primer_up,window.CHECK_n_inv_primer_down, window.ndeg_sites_up,window.ndeg_sites_down,window.n_variable_sites_in_window,window.primerup,window.primerdown]])+"\n")
			n_ok_windows+=1
with open(args.output_summary) as f:
	for window in windows:
		if not window.validation:
			output.write("\t".join([str(x) for x in [window.validation,window.start+1,window.end+1,window.CHECK_n_inv_primer_up,window.CHECK_n_inv_primer_down, window.ndeg_sites_up,window.ndeg_sites_down,window.n_variable_sites_in_window,window.primerup,window.primerdown]])+"\n")

output.close()

##Write new alignment
output=open(args.output_alignment,"w")
for key in dict_fasta.keys():
	output.write(">"+key+"\n"+dict_fasta[key]+"\n")
output.close()


#FINAL MESSAGE
print("\n\n\n\nSEARCH DONE:\nWe selected windows of %i bp  \nPrimers of %i bp with at least %i invariable sites at the end and max %i degenerated sites in each primer\n#windows passing criteria ( if any) at the beginning of the summary file\n" % (args.length_window,args.length_primer,args.n_inv_primer,args.n_deg_sites))
print("SUMMARY IN: "+args.output_summary+"\nNEW ALIGNMENT IN: "+args.output_alignment+"\nFOUND "+str(n_ok_windows)+" WINDOWS RESPECTING CRITERIA\nDONE")

