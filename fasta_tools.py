#!/usr/bin/env python
# Filename: fasta_tools.py
try:
	from Bio import SeqIO
	import windows_tools
except:
	raise Exception  ("could not import the modules required all or some of the functions : Bio or windows_tools")


###put the  the scaffolds in the assembly into a dictionnary 

file_scaf="/proj/b2010010/repos/assembly/fAlb15/linkage/fAlb15.chrom.strict.20140121.txt"
file_chrom="/proj/b2010010/repos/assembly/fAlb15/fAlb15.chrom.fa"
scaf_len="/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len"
class Assembly_scaff(object):
	"""class to store info on assembly from scaffolds information. needs a file where every line is a scaffold ordered along chromosome with three columns:
	chromosome, scaffold name, scaffold length  
	"""
	def __init__ (self,file_scaf):
		self.file_scaf=file_scaf
		self.dict_chrom={}
		self.dict_scaf={}
		with open (self.file_scaf,"r") as f:
			for line in f:
				self.dict_scaf[line.split()[1]]=windows_tools.WindowBed(seq=line.split()[1],start=0,length=line.split()[2])
				if line.split()[0] in self.dict_chrom.keys(): 
				 	self.dict_chrom[line.split()[0]].append(str(line.split()[1]))
				else:
				 	self.dict_chrom[line.split()[0]]=[line.split()[1]]
		self.dict_scaf_length={}
		with open(scaf_len) as f:
			for line in f:
				self.dict_scaf_length[line.split()[0]]=line.split()[1]

RNA_codon_table = {
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', 'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', 'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', 'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Urp',
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', 'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', 'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', 'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', 'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', 'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', 'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', 'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', 'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'}




#store info on assembly fAlb15
fAlb15=Assembly_scaff(file_scaf="/proj/b2010010/repos/assembly/fAlb15/linkage/fAlb15.chrom.all.20140121.txt")


dict_scaf_len={}
with open (scaf_len) as f :
	for line in f:
		dict_scaf_len[line.split()[0]]=line.split()[1]

scaffs=0
total_length=0
dict_chrom_len={}
for key in fAlb15.dict_chrom.keys():
	print key
	length_chrom=-5000# we add every scaf length_chrom + 5000... need to remove 5000 in the end cause the last scaffold is not followed by 5000 gaps
	for scaf in fAlb15.dict_chrom[key]:
		length_chrom += int(dict_scaf_len[scaf]) + 5000
		scaffs+=1
	dict_chrom_len[key] = length_chrom
	total_length+=length_chrom
#store scaf length for every scaffold in the assembly 

def parse_fasta_to_dict(filename,output_format="string"):
	''' parse a fasta file to a dictionnary of sequences with values being either seq objects from Biopython or string
	 and keys being sequence_name
	 ------------
	 Parameters:
	 filename="file"# name of the file you want to parse 
	 output_format="string|seq" #keys of  the dictionnary are either simple string or sequences objects from bipython'''
	fasta_sequences={}
	for seq_record in SeqIO.parse(filename, "fasta"):
	    if output_format=="string":
	    	fasta_sequences[seq_record.id]=str(seq_record.seq)
	    if output_format=="seq":
	    	fasta_sequences[seq_record.id]=seq_record
	print len(fasta_sequences),"sequences in the dictionnary"
	return fasta_sequences





def bed_from_fasta(fasta_file="file.fasta",output_file="file.bed"):
	"""create a bed file with 1 line per sequence in the fasta corresponding to the length of the sequence"""
	dict_seq = parse_fasta_to_dict(fasta_file,"string")
	output=open(output_file,"w")
	for key in dict_seq.keys():
		print key,len(dict_seq[key])
		output.write(key+"\t0"+"\t"+str(len(dict_seq[key]))+"\n") 
	output.close()

def exclude_lines_with_scaff(infile,outfile,list_scaff_to_exclude=fAlb15.dict_chrom["ChrZ"],nb_lines_log=10000):
	'''exclude lines with scaffolds specified in list_scaff_to_exclude list
	The scaff can be anywhere in the name'''
	outf=open(outfile,"wb")
	i=0
	with open(infile) as f:
		for line in f:
			i+=1
			if i%nb_lines_log == 0: print i," lines done"
			if not any([(to_exclude in line) for to_exclude in list_scaff_to_exclude]):
				outf.write(line)		
	outf.close()

