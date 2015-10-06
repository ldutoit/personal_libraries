#From ludovic Dutoit, July 14 2015
###Parameters

fasta = "masked_fasta.fa"
window_file = "windows.bed"
output_file = "genomic_features/GC/"




#Modules loading 
import os
os.system("module load  bioinfo-tools")
os.system("module load  biopython/1.56")

from Bio import SeqIO
from Bio.SeqUtils import GC






####FUNCTIONS



class WindowBed(object):
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
		if format_input!="BED": self.start=self.start-1# assume a one base inclusive and then it is a bed
		if self.end==None:
			self.end=int(self.start)+int(self.length)
		elif self.length==None:
				self.length=int(self.end)-int(self.start)
		if self.end<self.start: raise Exception("Error: Window end position is before start "+self.seq+" "+str(self.start)+" "+str(self.end))
		if self.length==0:  Warning("Error: Window has length 0 "+self.seq+" "+str(self.start)+" "+str(self.end))
		if  not strand in ["+","-",None]: raise Exception("Error: unauthorized strand: %s" % (str(self.strand)))
	def __repr__(self):
		return repr("    ".join([str(self.seq),str(self.start),str(self.end)]))

	def overlap(self,w2):
		"""for BED windows, return True if overlap
		>>> win=WindowBed(seq="Chr1",start=10,end=20)
		>>> win.overlap(WindowBed(seq="Chr1",start=10,end=20))
		True
		>>> win.overlap(WindowBed(seq="Chr1",start=20,end=30))
		False
		>>> win.overlap(WindowBed(seq="Chr2",start=10,end=30))
		False
		"""
		if self.seq != w2.seq: return False
		return  min(int(self.end), int(w2.end)) - max(int(self.start), int(w2.start)) > 0

	def overlap_length(self,w2):
		"""for BED windows, return length of overlap
		>>> win=WindowBed(seq="Chr1",start=10,end=20)
		>>> win.overlap_length(WindowBed(seq="Chr1",start=10,end=20))
		10
		>>> win.overlap_length(WindowBed(seq="Chr1",start=20,end=30))
		0
		>>> win.overlap_length(WindowBed(seq="Chr2",start=10,end=30))
		0
		"""
		if self.seq != w2.seq: return 0
		return max(0, min(self.end, w2.end) - max(self.start, w2.start))
	
	# def overlap_listwindows(self,list_windows,out="nb",verbose=True):# retun number of hits or number of bp
	# 	"""for BED windows, return length of overlap or nb hits to a list of windows
	# 	>>> win=WindowBed(seq="Chr1",start=10,end=20)
	# 	>>> list_win=[WindowBed(seq="Chr1",start=10,end=20),WindowBed(seq="Chr2",start=10,end=20)]
	# 	>>> win.overlap_length(WindowBed(seq="Chr1",start=10,end=20))
	# 	10
	# 	>>> win.overlap_listwindows(list_win,verbose=False)
	# 	10
	# 	"""
	# 	win=self
	# 	print win
	# 	nhits=0
	# 	bp=0
	# 	for win_to_check in list_windows:
	# 		ov_length=self.overlap_length(win_to_check)
	# 		if ov_length!=0:
	# 			nhits+=1
	# 			bp+=ov_length
	# 			if verbose==True: print "OVERLAP: ", self,win_to_check
	# 	if out=="nb":
	# 		return nhits
	# 	elif out=="bp":
	# 		return bp
	# 	else:
	# 		raise Exception("invalid out arguments (choose 'bp'|'nb')")

	def overlap_listwindows(self,list_windows,out="nb",verbose=True, col_interest = 4, col_weight  = None ):# retun number of hits or number of bp
		"""for BED windows, return length of overlap or nb hits to a list of window.  out can be bp or nb (number of overlaps or of bpair overlap)
		NOTE: out ="weighted_average"-It can also return average value for a given column (col_interest = INT, 1 based) in the list of windows overlapping with the original window. It is weighted by the nb of bp overlapping for every overlap.
		One can give another column as the weight using col_weight , if None.. the function use the number of bp overlapping to weight
		It would then the average value as well as the nb of bp overlap
>>> win=WindowBed(seq="Chr1",start=10,end=20)
>>> list_win=[WindowBed(seq="Chr1",start=10,end=20),WindowBed(seq="Chr2",start=10,end=20)]
>>> win.overlap_length(WindowBed(seq="Chr1",start=10,end=20))
10
>>> win.overlap_listwindows(list_win,verbose=False)
'Chr1    10    20'
1
>>> win.overlap_listwindows(list_win,verbose=False,out="bp")
'Chr1    10    20'
10
    	"""
		win=self
		print win
		nhits=0
		bp=0
		if out == "weighted_average": weight = unweighted_sum = 0 # for weighted average
		for win_to_check in list_windows:
			ov_length=self.overlap_length(win_to_check)
			if ov_length!=0:
				nhits+=1
				bp+=ov_length
				if out =="weighted_average":
					value = float(win_to_check.allcols[col_interest-4])
					if col_weight != None : 
						temp_weight = float(win_to_check.allcols[col_weight-4])
						weight+=temp_weight
					else:
						temp_weight = ov_length
					unweighted_sum += value * int(temp_weight)
				if verbose==True: print "OVERLAP: ", self,win_to_check
		if out=="nb":
			return nhits
		elif out=="bp":
			return bp
		elif out=="weighted_average":
			if bp != 0:
				return [unweighted_sum/weight,bp]
			else:
				return ["NA",0]
		else:
			raise Exception("invalid out arguments (choose 'bp'|'nb' | 'weighted_average')")

	def write(self,filehandle,optional_list=[],allcols=True):
		""" write a window to a filehandle,optional_list ist a list that will be appended to the window, allcols decide if all win.allcols is printed too as expanded columns
		>>> win=WindowBed(seq="Chr1",start=10,end=20)
		>>> output = open ("temp","w")
		>>> win.write(output,optional_list =["test"],allcols = True) # should just contain test as there is no allcols
		>>> win.write(output,allcols = False) ### just window
		>>> win.allcols = ["col1",2] # extracols who are not necessarily strting
		>>> win.write(output)
		>>> output.close()
		>>> filecmp.cmp("temp","test_files/windows_write.txt") #check that these test generate the same file that a precompiled file
		True
		>>> os.remove("temp")
		"""
		standard=[str(info) for info in [self.seq,self.start,self.end]]
		if allcols==True and "allcols" in self.__dict__ : others = [str(info) for info in self.allcols] + optional_list
		else : others= optional_list
		all_print = "\t".join(standard+others)
		filehandle.write(all_print+"\n")

# # 	# def get_seq(self,fasta):
# # 	# 	"""Extract the DNA sequence corresponding to your window in a file
# # 	# 	>>> win=WindowBed(seq="Chr1",start=10,end=20)
# # 	# 	>>> fasta_sample=
# # 	# 	"""
# # 	# 	sequences=fasta_tools.parse_fasta_to_dict(fasta,output_format="seq")
# # 	# 	self.sequence=sequences[self.seq].seq[self.start:self.end]




class Bed(object):
	def __init__ (self,path,seqcol,startcol,endcol=None,lengthcol=None,strandcol=None,format_input="BED",allcols=True,skip_nlines=0):
		"""
		a class for Bedfile
		col position are 0 based. Will create a list of windows bed based at Bed.windows
		If format_input!=BED, it assumes a 1 based inclusive format and transfor it into bed
		allcols=TRUE|FALSE is a list containing all the columns of the file that are not stored in the rest
		skip_nlines number of header lines to skip, (inclyude line starting by # that are skipped anyway. 
		
		>>> a = Bed("test_files/windows_write.txt",0,1,2)
		>>> a
		'Bed from test_files/windows_write.txt contains: 3 windows'
		>>> a.windows
		['Chr1    10    20', 'Chr1    10    20', 'Chr1    10    20']
		>>> b = Bed("test_files/windows_write.txt",seqcol=0,startcol=1,endcol=2,format_input="noBed") 
		>>> b.windows
		['Chr1    9    20', 'Chr1    9    20', 'Chr1    9    20']
		"""
		self.path=path
		self.format=format
		#add the values
		if  all(v is None for v in [endcol,lengthcol]): 
			raise Exception("Window definition error: end and length have to be defined, provide at least one of them")
		if format!="BED": "python WindowBed will be zero based. Assume that the input format is vcf/GFF like , 1 based and inclusive"
		windows=[]
		with open(path) as f:
			i=0
			for line in f:
				i+=1
				if not line.startswith("#"):
					if i>skip_nlines:
						info=line.split()
						if endcol!=None:
							win=WindowBed(seq=info[seqcol],start=info[startcol],end=int(info[endcol]),format_input=format_input)
						elif lengthcol!=None:
							win=WindowBed(seq=info[seqcol],start=info[startcol],length=int(info[lengthcol]),format_input=format_input)
						if strandcol!=None:
							win.strand=info[strandcol]
						if allcols==True:
							toadd=[j for i,j in enumerate(info) if not i in [seqcol,startcol,endcol,lengthcol]] # add only the non primary argument
							win.allcols=toadd
						windows.append(win)
		self.windows=windows
	def __repr__(self):
		return repr("Bed from "+self.path+" contains: "+str(len(self.windows))+" windows")


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


#*CORE PROGRAM
# read fasta
sequences_masked=fasta_tools.parse_fasta_to_dict(fasta,output_format="seq")

windows = wt.Bed(window_file,0,1,endcol=2).windows
output=open(output_file,"w")
for win in windows:
	win.fa = str(sequences_masked[win.seq][win.start:win.end].seq)
	countnoN = win.length-win.fa.count("N") # windows length - N
	countGC =win.fa.count("G") + win.fa.count("C")
	win.GC = countGC/float(countnoN)
	print win.seq,win.start,win.end,win.GC
	output.write("\t".join( [win.seq,str(win.start),str(win.end),str(win.GC)])+"\n")
output.close()


