
###Store all windows
import windows_tools as wt
import os


originalCDS = "beds/Ficedula_albicollis.fAlb15.e73_clean_sortedCDSonly.bed"
genewindows = "beds/Ficedula_albicollis.fAlb15.gene_windows.bed"

CDSwithgeneinfo = wt.Bed(originalCDS,0,1,2)
genewindows =wt.Bed(gene_bed_file,0,1,2)
if not os.path.exists("tempconcat"):
	os.system("cat beds/CDSbeds/* > tempconcat")

bedobject = wt.Bed("tempconcat",0,1,2)
dict_per_scaf = {}
for win in bedobject.windows:
	if not win.seq in dict_per_scaf.keys():
		dict_per_scaf[win.seq] = [win]
	else :
		dict_per_scaf[win.seq].append(win)



###test for overlapping windows
overlaps = []# put pairs in list, remove genes
for key in dict_per_scaf.keys():
	print key
	listwin = dict_per_scaf[key]
	while len(listwin)>0:
		overlappingwindows = listwin[0].overlap_listwindows(listwin,out="list",verbose=False)
		if len(overlappingwindows)>1:
			overlaps.append([listwin[0],overlappingwindows[0]])
		listwin.remove(listwin[0])



#Do the remvoing part
geneset = set()
for item in genewindows.windows:
	a = item.overlap_listwindows(genewindows.windows,out="list")
	if len(a)>1: 
		for overlap in a:
			geneset.append(a)
weird2 = []
geneset2 = set()
i=0
for item in overlaps:
	i+=1
	print i
	a = item[0].overlap_listwindows(genewindows.windows,out="list")
	if len(a)>0: 
		geneset2.add(a[0].allcols[0])
	b = item[1].overlap_listwindows(genewindows.windows,out="list")
	if len(b)>0: 
		geneset2.add(b[0].allcols[0])
	if len(a) ==0:
		weird2.append(item[0])
	if len(b) ==0:
		weird2.append(item[1])
		