#!/usr/bin/env python 2
import random,vcf
import vcf_tools as vcf_t
import windows_tools
import sys
import gc

# run bootstraps pi bootraps just run pi bootstrap several time, it is highly dependant on the funciton pi_double_vcf_fix_nsites

def pi_bootstrap(windows,vcf_file_snps,vcf_file_allsites,mincov=5,maxcov=100,nbMark = 10,nbInd=2,minNsites = 400,maxNsites = 800,bgzip=True) :
	#pick random nbInd individuals
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file_snps, compressed=bgzip, prepend_chr="False", strict_whitespace=False)
	inds = random.sample(input_vcf.samples,nbInd)
	#sample marker
	markers=[]
	sum_pairwise=0
	sum_nsites=0
	while len(markers) < nbMark:
		print "\nnb_markers:"+str(len(markers))
		marker = random.choice(windows)#pick a marker
		print marker,marker.allcols
		if not marker in markers: 
			counts = vcf_t.pi_double_vcf_fix_nsites(vcf_file_snps,vcf_file_allsites,marker.seq,marker.start,marker.end,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=True,min_nsites=minNsites,max_nsites=maxNsites,called=True)
			print counts
			print markers
			if counts[0]!= "NA":
				sum_pairwise +=float(counts[0]) * float(counts[1]) #pi*nsites_ok
				sum_nsites += float(counts[1])#nsites_ok  
				markers.append(marker)
				print sys.getsizeof(markers)
	return sum_pairwise/sum_nsites

def pi_bootstraps(bed_file,vcf_file_snps,vcf_file_allsites,nb_bootstraps=1000,seqcol=0,startcol=1,endcol=2,output_file="test_bootstraps.txt",
		mincov=5,maxcov=100,nbMark = 10,nbInd=2,minNsites = 400,maxNsites = 800,bgzip=True):
	bootstraps=[]
	windows = windows_tools.Bed(bed_file,0,1,2).windows
	print "markers loaded"
	print sys.getsizeof(windows)
	while len(bootstraps) < nb_bootstraps:
		pi=pi_bootstrap(windows,vcf_file_snps,vcf_file_allsites,mincov=mincov,maxcov=maxcov,nbMark = nbMark ,nbInd= nbInd ,minNsites = minNsites,maxNsites = maxNsites,bgzip=bgzip)
		bootstraps.append(pi)
		print "bootstrap nr",len(bootstraps),pi
		print sys.getsizeof(pi)
		print sys.getsizeof(bootstraps)
		gc.collect()
	#output
	output=open(output_file,"w")
	for bootstrap_pi in bootstraps:
		output.write(str(bed_file)+"\t"+str(nbMark)+"\t"+str(nbInd)+"\t"+str(maxNsites)+"\t"+str(bootstrap_pi)+"\n")
	output.close()


