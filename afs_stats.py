 #!/usr/bin/env python 2
#Filename: afs_stats.py

import numpy 

def asf_basic_stats(afs):
	''' 
	Calculate Tajima's D,theta watterson and nucleotide diversity from the allele frequency spectrum.
	Input: the site frequency spectrum as a list for x chromosmes:
	[1/x, 2/x, 3/x....x-1/x, total valid sites]
	It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list. 
	It returns, pi, theta watterson and TajD
	Example:
	>>> asf_basic_stats([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[1]: 0.093333333333333338
	>>> asf_basic_stats([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[2]: 0.093333333333333338
	'''
	mythetaW = afs_thetaW(afs) # calculate Theta Watterson
	mypi = afs_nucl_diversity(afs) # calculate nucleotide diversity
	varD =  variance_D(afs) # calculate the variance for d
	TajD = (mypi-mythetaW) / varD # calculate Tajima's D
	print " returning pi, theta and TajD...."
	return mypi, mythetaW, TajD


def afs_thetaW(afs):
	''' 
	Calculate theta watterson from the allele frequency spectrum.
	Input: the site frequency spectrum as a list for x chromosmes:
	[1/x, 2/x, 3/x....x-1/x, total valid sites]
	It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list.
	Example:
	>>> afs_thetaW([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[1]: 0.08759124087591241
	>>> afs_thetaW([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[2]: 0.08759124087591241
	'''
	nchrom = len(afs) # the number of chromomosomes
	seg_site = sum([int(x) for x in afs[:-1]]) # number of segregating sites
	n_total = float(afs[-1]) # total number of sites
	harmonic_mean = sum([1/float(x) for x in range(1,nchrom)]) # harmonic mean 
	theta_watterson =(seg_site/harmonic_mean) / n_total
	print "assuming",nchrom,"chromosomes",seg_site,"variable sites"
	return theta_watterson


def afs_nucl_diversity(afs):
	''' 
	Calculate nucleotide diveristy (pi) from the allele frequency spectrum.
	Input: the site frequency spectrum as a list for x chromosmes:
	[1/x, 2/x, 3/x....x-1/x, total valid sites]
	It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list.
	Example:
	>>> afs_nucl_diversity([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[1]: 0.093333333333333338
	>>> afs_nucl_diversity([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
	Out[2]: 0.093333333333333338
	'''
	n = len(afs) # number of chromosomes
	sf = numpy.array([float(x) for x in afs[:-1]]) # allelc frequency spectrum
	tot = afs[-1] # number of valid sites
	pi = sf*2* numpy.array(range(1,n)) * numpy.array(range(n-1,0,-1)) 
	pi2 = (sum(pi /(n*(n-1)))) /float(tot)
	return pi2


def variance_D(afs):
	'''Calculate Tajima's D,theta watterson and nucleotide diversity from the allele frequency spectrum.
	Input: the site frequency spectrum as a list for x chromosmes:
	[1/x, 2/x, 3/x....x-1/x, total valid sites]
	It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list. 
	'''
	n = len(afs) # number of chromosomes
	a1 = sum([1.0/x for x in range(1,n)])
	b1 = (n+1.)/(3*(n-1)) 
	c1 = b1-( 1/a1)
	e1 = c1/a1

	a2 = sum([1.0/(x**2.) for x in range(1,n)])
	b2 = ( 2*(n**2.+n+3.)  )/ ( (9.0*n) *  (n-1.)      )
	c2 = b2 - ((n+2.)/(a1*n)) + (a2)/(a1**2.)
	e2 = c2 /((a1**2.)+(a2))
	
	s = sum([int(x) for x in afs[:-1]]) # number of segregating sites
	varD = (((e1*s)+(e2*s*(s-1)))**0.5)
	return varD

