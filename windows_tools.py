#!/usr/bin/env python2
# Filename: windows_tools.py
#I try to import the modules several time because it appeared several times that the import fail due to other parralel import when running several parralel instances of the same script.
for i in range(0,1000):
    while True:
        try:
			import os,sys,vcf,pysam
			import numpy as np 
			import filecmp 
        except :
            raise Exception  ("could not import the modules required all or some of the functions : pysam\nnumpy\nos\nvcf\nrandom"+str(i))
        break
        
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
	def overlap_window(self,w2):
		"""for BED windows, return a window with the coordinates of the overlap

		"""
		assert self.overlap(w2)# they have to overlap
		return WindowBed(self.seq, max(self.start, w2.start), min(self.end, w2.end)) # from the biggest start to the smaller end
	def gap (self,w2,out ="nb"):
		""" 
		For BED windows, return the gap gap between two windows.
		out = "nb"| "window" # return the nb of the gap b default but can also return coordinates of the gap
		>>> win2 = WindowBed("a",100,200)
		>>> win1 =  WindowBed("a",250,300)
		>>> win1.gap(win2)
		50
		"""
		assert self.seq  == w2.seq  # same sequence for both windows
		assert not self.overlap(w2)  # there is a gap between them
		if out == "nb":
			return  max(self.start,w2.start) - min(self.end,w2.end)
		elif out == "window":
			return WindowBed(self.seq,min(self.end,w2.end),max(self.start,w2.start))
		else:
			raise Exception("invalid out argument, should be 'nb' or 'window'")

	def no_overlap_listwindows(self,list_windows,verbose=True ):# 
		''' check that the window has no overlap with a list of windows, return True if no overlap or False otherwise'''
		win=self
		if verbose==True: print( win)
		nhits=0
		bp=0
		hits = []
		for win_to_check in list_windows:
			if win.overlap(win_to_check): return False
		return True
	def overlap_listwindows(self,list_windows,out="nb",verbose=True, col_interest = 4, col_weight  = None ):# retun number of hits or number of bp
		"""for BED windows,  tools to overlap a window a list of windows. SEE ovarlap_two_list_windows for directly compared two windows list
			
		list_windows  # the list of windows you want to check for overlaps 

		out can be bp or nb (number of overlaps or of bpair overlap)
		NOTE: out ="weighted_average"-It can also return average value for a given column (col_interest = INT, 1 based) in the list of windows overlapping with the original window. It is weighted by the nb of bp overlapping for every overlap.
		One can give another column as the weight using col_weight , if None.. the function use the number of bp overlapping to weight
		It would then the average value as well as the nb of bp overlap
		out = nb,bp, weighted_average or list

		>>> win=WindowBed(seq="Chr1",start=10,end=20)
		>>> list_win=[WindowBed(seq="Chr1",start=10,end=20),WindowBed(seq="Chr2",start=10,end=20)]
		>>> win.overlap_length(WindowBed(seq="Chr1",start=10,end=20))
		10
		>>> win.overlap_listwindows(list_win,verbose=False)
		1
		>>> win.overlap_listwindows(list_win,verbose=False,out="list")
		['Chr1    10    20']
		>>> list_win[0].allcols = [10,10]###to test weighted average
		>>> list_win.append(WindowBed(seq="Chr1",start=5,end=15))
		>>> list_win[2].allcols = [5,10]
		>>> list_win[1].allcols = [10,20]
		>>> win.allcols = [1,2]
		>>> win.overlap_listwindows(list_win,verbose=False,out="weighted_average",col_interest = 4)
		[8.333333333333334, 15]
		>>> win.overlap_listwindows(list_win,verbose=False,out="weighted_average",col_interest = 4,col_weight  = 5)
		[7.5, 15]
    	"""
		win=self
		if verbose==True: print( win)
		nhits=0
		bp=0
		hits = []
		if out == "weighted_average": weight = unweighted_sum = 0 # for weighted average
		for win_to_check in list_windows:
			ov_length=self.overlap_length(win_to_check)
			if ov_length!=0:
				hits.append(win_to_check) 
				nhits+=1
				bp+=ov_length
				if out =="weighted_average":
					if verbose==True: print( win_to_check)
					value = float(win_to_check.allcols[col_interest-4])
					if col_weight != None : 
						temp_weight = float(win_to_check.allcols[col_weight-4])
					else:
						temp_weight = ov_length
					weight+=temp_weight
					unweighted_sum += value * float(temp_weight)
				if verbose==True: print ("OVERLAP: ", self,win_to_check)
		if out=="nb":
			return nhits
		elif out=="bp":
			return bp
		elif out=="weighted_average":
			if bp != 0:
				if verbose==True: print (" [weighted average ,", "total_weight]")
				return [unweighted_sum/weight,bp]
			else:
				return ["NA",0]
		elif out=="list":
			return hits
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

	# def write(self,path="windows.bed",format="BED",extrainfo="all",headers=False,col3="end",strand=False):
	# 	"""Method to write quickly a window file from a 'Bed object. 
		
	# 	path="windows.bed"# path to which the window file is outputted
		
	# 	format="BED"# If not BED, assume a 1 based inclusive format GFF type
		
	# 	extrainfo="all"|None|[poscol,poscol]#to rite all the extra info stored per window. All write all extra columns, None write none, 
	# 	and otherwise one has to use a list of indices (0 based)
		
	# 	headers=False#A list that will serve as header for the file a # will be added at the start
	# 	col3="end"|"len"# define if the thirdcol is the length or the end of the interval. end by default
		
	# 	strand=True output the strand as the 4th colum

	# 	"""
	# 	if format!="BED": print "assume that python as zero based windows and print a 1 based format inclusive(GFF type)  "
	# 	#opening thw output
	# 	output=open(path,"w")
	# 	if headers!=False: output.write(headers+"\n")#writing the potential header
	# 	if extrainfo=="all": 
	# 		extrainfo=range(0,len(self.windows[0].allcols))
	# 	elif extrainfo is not list and extrainfo is not "all":
	# 		raise Exception("extrainfo='all'|None|[poscol,poscol]#to rite all the extra info stored per window. All write all extra columns, None write none, \
	# 			and otherwise one has to use a list of indices (0 based)")
	# 	for window in self.windows:
	# 		if format!="BED": window.start+=1
	# 		#print the column order
	# 		if col3=="len":
	# 			basics=([str(window.seq),str(window.start),str(window.length)])
	# 		elif col3=="end":
	# 			basics=([str(window.seq),str(window.start),str(window.end)])
	# 		else:
	# 			raise Exception("col3 argument  must be either'len' either 'end'")
	# 		if strand==True:
	# 			basics+=window.strand
	# 		extrainfo_to_write=[]
	# 		if extrainfo!=None:
	# 			for colnumber in extrainfo:
	# 				if len(window.allcols[colnumber])!=0:#print the annotation to the file
	# 					extrainfo_to_write+=[str(window.allcols[colnumber])]
	# 		list_towrite=basics+extrainfo_to_write+["\n"]
	# 		output.write("\t".join(list_towrite))
	# 	output.close()
	# 	print "wrote ",len(self.windows)," to ",path

def scaf_length_to_windows_bed(scaflength_file="/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len", window_size=10000,output_file="windows.bed",cover_all_scaff=True):
	'''A function that return a bed file of discret windows along scaffolds
	Parameters:
	scaflength_file="/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len" a file containing two columns, col 1 scaffnames, col 2 length
	window_size=100 in bp
	output_file="windows.bed" the name of the output file containing windows
	cover_all_scaff=True if True the program make a shorter ewindow to finish the scaffold even if it is shorter, if False only windows of  exact window size
	>>> scaf_length_to_windows_bed(scaflength_file="/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len", window_size=10000,output_file="temp.bed",cover_all_scaff=True)
	Some windows will be smaller than the desired size as we cover all scaffolds (cover_all_scaff == True)
	>>> sum(1 for line in open('temp.bed'))
	131139
	>>> scaf_length_to_windows_bed(scaflength_file="/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len", window_size=10000,output_file="temp.bed",cover_all_scaff=False) # a file w
	>>> sum(1 for line in open('temp.bed'))
	109296
	'''
	## extract a list of scaffolds and scaffold length
	if cover_all_scaff==True: print "Some windows will be smaller than the desired size as we cover all scaffolds (cover_all_scaff == True)"
	scaffolds=[]
	if not os.path.isfile("/proj/b2010010/private/assembly/nobackup/ScaffLengths/fAlb15.len") : raise Exception ("scaffold file:" ,scaflength_file, "does not exist")
	with open(scaflength_file) as f:
		for line in f:
			info=line.split()
			scaffolds.append(info)
	# go scaffold by scaffold and output windows to the bed file
	output=open(output_file,"w")# open output file
	for scaff,scafflength in scaffolds:# # for every scaffold, scaffold length
		#print scaff 
		pos=0 #allow to start at 0 
		while pos<int(scafflength)-int(window_size):# as long as we can make window of window size in the scaffold
			start=pos# start position (first window start at 0)
			end=pos+window_size
			pos=end
			output.write("%s\t%i\t%i\n" % (scaff,start,end))
		if pos < int(scafflength) and cover_all_scaff==True:
			#print "last window"
			start=pos
			end=int(scafflength)+1
			output.write("%s\t%i\t%i\n" % (scaff,start,end))
	output.close()


def intersect_per_windows(window_file,annotation_file,output_folder):
	''' take a window file (bed) and intersect it with an another window file (GFF/GTF/BED). Report one file per window in the window file conataining all the regions that intersect with the window_file.
	NOTE automatically correct for GTF strting 1
	Parameters:
	window_file# a windows file for which you want a file per intersect
	annotation_file# a file for which you want to report the entries per window of window file
	output_folder#the folder in which you will have a folder per window in the window file
	'''
	# gff, gtf or bed?
	file_extension=annotation_file.split(".")[len(annotation_file.split("."))-1]
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)
	else:
		raise Exception("The output folder already exists!")
	if file_extension in ["gtf","gff","bed"]: # create a fake bed file
		os.system("intersectBed -a "+window_file  +" -b "+ annotation_file +" -wa -wb  > "+window_file+"temp.bed")
		with open(window_file) as f:
			for line in f:
				os.system("cat "+window_file + "temp.bed | grep -e '^"+line.split()[0]+"\s*"+line.split()[1]+"\s*"+line.split()[2]+"' | sed -e  's/^"+line.split()[0]+"\s*"+line.split()[1]+"\s*"+line.split()[2]+"\s*//g' > "+output_folder+"/"+"_".join(line.split())+".bed" )
	else:
		"Annotations file can only handle a .gff , .gtf or .bed"
	os.remove(window_file + "temp.bed")
	
def intersect_per_window_vcfgz(vcfzipped_file,window_file,output_folder,rm=False):
	'''similar as intersect_per_windows but output a vcf file per line int he window_filke file. You  use tyo only output rmasked region by setting to rm to a bedffile of repeat. 
	It should be a bedfile of repeats and not a bedfile of regions non-repeatedcan
	NOTE automatically correct for GTF strting 1
	Parameters:
	window_file# a windows file for which you want a file per intersect
	window_file# a file for which you want to report the entries per window of window file
	output_folder#the folder in which you will have a folder per window in the window file
	rm# False, if stated otherwise will i
	'''
	os.mkdir(output_folder)
	with open(window_file) as f:
		for item in f:
			print "Adding a repeat mask"
			#print item.split()
			info=item.split()
			#print ("tabix "+VCF_allsite_file+" "+info[0]+":"+str(int(info[1])+1)+"-"+str(int(info[2]))  +" >  "+folder_vcf_files_all_sites_per_markers+"/"+info[0]+"_"+str(int(info[1])+1)+"_"+str(int(info[2]))+"_"+str(info[3])+".vcf")
			if rm==False:
				os.system("tabix "+vcfzipped_file+" "+info[0]+":"+str(int(info[1])+1)+"-"+str(int(info[2]))  +" >  "+output_folder+"/"+info[0]+"_"+str(int(info[1])+1)+"_"+str(int(info[2]))+".vcf")
			else:
				os.system("tabix "+vcfzipped_file+" "+info[0]+":"+str(int(info[1])+1)+"-"+str(int(info[2]))  +" | intersectBed -a - -b "+rm+  " -v  >  "+output_folder+"/"+info[0]+"_"+str(int(info[1])+1)+"_"+str(int(info[2]))+".vcf")

def gtftobed(infile,outfile,keep_annotations=False):
	'''transform a gtf or a gtf to a bed file ut loose annotation'''
	output=open(outfile,"w")
	with open(infile) as f:
		for line in f:
			info=line.split("\t")
			if keep_annotations==False:
				output.write(str(info[0])+"\t"+str(int(info[3])-1)+"\t"+str(info[4])+"\n")
			else:
				output.write(str(info[0])+"\t"+str(int(info[3])-1)+"\t"+"\t".join(info[4:len(info)])+"\n")
	output.close()


def sortabed(infile):
	''' sort a bed file and output the same infile than outfile'''
	os.system("sort -k1,1 -k2,2g -o temp.bed " +infile)
	os.remove(infile)
	os.rename("temp.bed",infile)

def subtractBed(infile,outfile,regionstosubtract):
	''' simplifiedsubtractBed from bedtools'''
	os.system("subtractBed -a "+ infile + " -b " + regionstosubtract + " | sort -k1,1 -k2,2g - >"  + outfile)


def subtractBed_recursively(folder_infiles=".",regionstosubtract="",outfolder="subtractedBeds",prefixoutfiles="subtracted."):
	''' apply recursively  simplified subtractBed from bedtools to allk the bedtools in a folder. Direct output into a new folder'''
	os.mkdir(outfolder)
	for item in os.listdir(folder_infiles):
		if item.endswith(".bed"):
			print str(item)
			subtractBed(infile=item,outfile=str(outfolder)+"/"+prefixoutfiles+item,regionstosubtract=regionstosubtract)
		else:
			print str(item) + " is not a '.bed' file... skipping"


def get_headers(file):
	''' get the header of a VCF file and return it in a list'''
	with open(file) as f:
		for line in f:
			if not line.startswith("##"):
				if line.startswith("#"):
					headers=line.split()
					return headers
					break


def splitBed(bedfile,nb=10,per="bp",output_folder="splitted",prefix="splitted_"):
	'''split a bedfile into a new folder into "nb" files in alphabetical order (up to 10000) per "bp" or "num_win"
	'''
	assert per=="bp" or per=="nb_win","invalid per argument, should be bp or num_win"
	windows = Bed(bedfile,0,1,endcol=2).windows
	os.mkdir(output_folder)
	file_num=1
	buffer0="".join((str(i) for i in np.repeat(0,5-len(str(file_num)))))
	output=open(output_folder+"/"+prefix+buffer0+str(file_num)+".bed","w")
	if per=="bp":
		bp=0
		for win in windows:
			bp+=win.length
			win.write(output)
			if bp>=nb:
				output.close()
				file_num+=1
				buffer0="".join((str(i) for i in np.repeat(0,5-len(str(file_num)))))
				output=open(output_folder+"/"+prefix+buffer0+str(file_num)+".bed","w")
				bp=0
	if per=="nb_win":
		nb_win=0
		for win in windows:
			win.write(output)
			nb_win+=1
			if nb_win>=nb:
				output.close()
				file_num+=1
				buffer0="".join((str(i) for i in np.repeat(0,5-len(str(file_num)))))
				output=open(output_folder+"/"+prefix+buffer0+str(file_num)+".bed","w")
				nb_win=0
	output.close()
	print file_num," files in", output_folder

def bed_calculate_densities(bed,intersectwith,output_file):
	'''calculate the densities of each window in bed for regions in bed_denstities
	and append it as a last column in output_file
	NOTE Careful with redundancy, crete bugs!'''
	windows = Bed(bed,0,1,endcol=2).windows
	tointersectwith =  Bed(intersectwith,0,1,endcol=2).windows
	output=open(output_file,"w")
	i=0
	for win in windows:
		i+=1
		print win,i
		total = win.overlap_listwindows(tointersectwith,out="bp",verbose=False)
		win.dens =  str(total/float(win.length))
		win.write(output,[win.dens])
	output.close()

def bed_calculate_weighted_average(bed,intersectwith,output_file,col_interest=4,col_weight= None):
	'''calculate the weighted average per basepair overlapping of each window in bed for col_intersect(1 based)  in intersectwith
	and append two columns in output file, average value and bp overlapping
	NOTE Careful with redundancy, create bugs!
	Can also use a col_weight (1 based) instead of a window list'''
	windows = Bed(bed,0,1,endcol=2).windows
	tointersectwith =  Bed(intersectwith,0,1,endcol=2).windows
	output=open(output_file,"w")
	i=0
	for win in windows:
		i+=1
		print win,i
		win.weighted_average = win.overlap_listwindows(tointersectwith,out="weighted_average",verbose=False,col_interest=col_interest,col_weight  = col_weight )
		win.write(output,[str(win.weighted_average[0]),str(win.weighted_average[1])])
	output.close()

import doctest
doctest.testmod()


def overlap_two_list_windows(list1,list2):
	'''return a list of pairs overlaps between two list of WindowBed objects [[win1,win2],[win1,win2]]'''
	# create a dictionnary of list2 by sequence
	window_dict= {}
	for win in list2:
		if win.seq in window_dict.keys():
			window_dict[win.seq].append(win)
		else:
			window_dict[win.seq] = [win]
	#create a list of overlap using the windows
	overlaps = []
	i=0
	for win in list1:
		i+=1
		print "win number",i
		a = win.overlap_listwindows(window_dict[win.seq],out="list")
		if len(a)>1: 
			for overlap in a:
				overlaps.append([win,overlap])
	return overlaps

