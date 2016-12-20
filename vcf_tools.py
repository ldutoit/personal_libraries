#!/usr/bin/env python 2
#Filename: vcf_tools.py
try :
	import pysam,numpy,os,vcf,random,sh,gzip
except:
	raise Exception  ("could not import the modules required all or some of the functions : pysam\nnumpy\nos\nvcf\nrandom\nsh")

#sample_vcfgz=os.path.dirname(__file__)+"/test_files/sample_allsites.vcf.gz" # alternatively sample_vcfgz= "/home/ludovic/repos/personal_libs/test_files/sample_allsites.vcf.gz"
#sample_vcf=os.path.dirname(__file__)+"/test_files/sample_allsites.vcf.gz" # alternatively sample_vcf = "/home/ludovic/repos/personal_libs/test_files/sample_allsites.vcf"



def extract_pi_double_vcf_bed(vcf_snps,vcf_allsites,bed,output,mincov=5,maxcov=10000,inds="all",bgzip=True,min_nsites=10000,called=True,nb_ind_with_min_cov="all",uppmax=False):
	'''append pi information window by window to a bed file, for parameters information refers to the functions called

	vcf_snps # the path the the snps files
	vcf_allsites # the path to the allsites files
	bed # the bed containing the window for which you want the sites
	output ##the path to an output file
	mincov=5 #the minimum coverage for a given site that all individuals must have for the site to be considered but see  nb_ind_with_min_cov )
	maxcov=10000 #the maximum coverage, if any individual in "inds" have more, the site is excluded
	inds="all" # a list of individuals you want to estimate nucleotide diversity for, by default all the individuals in the dataset
	bgzip=True #wether the vcffiles are gzipped or not
	min_nsites=10000 # the minimum number of valid sites in any given window for pi to be estimated, if less return NA
	called=True # all individuals should be assigned a genotype (for SNPs)
	nb_ind_with_min_cov="all"  # the number of individuals in "inds" that should have the mincov for a given site, by default all individuals#
	'''
	if uppmax == True:
		scratch_folder = os.environ.get("SNIC_TMP")
		print "copying the vcfs to the scratch folder "+scratch_folder
		###copy thr vcf
		if not os.path.exists(scratch_folder+"/"+os.path.basename(vcf_snps)):
			os.system("cp "+ vcf_snps+" $SNIC_TMP/" )
			print "vcf_snps copied"
		if not os.path.exists(scratch_folder+"/"+os.path.basename(vcf_allsites)):
			os.system("cp "+ vcf_allsites+" $SNIC_TMP/" )
		print "copy complete"
		###copy the index files
		if not os.path.exists(scratch_folder+"/"+os.path.basename(vcf_allsites+".tbi")):
			os.system("cp "+ vcf_allsites+".tbi $SNIC_TMP/" )
		if not os.path.exists(scratch_folder+"/"+os.path.basename(vcf_snps+".tbi")):
			os.system("cp "+ vcf_snps+".tbi $SNIC_TMP/" )
		vcf_snps = scratch_folder+"/"+os.path.basename(vcf_snps)
		vcf_allsites = scratch_folder+"/"+os.path.basename(vcf_allsites)
	output=open(output,"w")
	with open(bed) as f:
		for line in f:
			info = line.split()
			count = pi_double_vcf(vcf_snps,vcf_allsites,line.split()[0],int(line.split()[1])+1,int(line.split()[2]),mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,min_nsites=min_nsites,called=True,nb_ind_with_min_cov=nb_ind_with_min_cov)
			info.append(str(count[0]))#sum_pairwise
			info.append(str(count[1]))#number of ok sites
			info.append(str(count[3]))#number of snps
			print "\t".join(info)+"\n"
			output.write("\t".join(info)+"\n")
	output.close()



def extract_features(vcf_file,feature,chrom,start,end,sampling_frequency=1,missing_data="None",bgzip=True,write=True,outfile=""):
	"""extract sample features from VCF format (i.e. DP/GT etc) on a table : Chrom,Pos,ind1...indN
	
	vcf_file # the vcf file from which to extract the data
	feature # DP:GT or any other abbreviation present in some lines of t
	he VCF. all ind will have missing data value for lines where the feature is not present
	chrom #OPTIONAL specify a sequence for which to extract the featurE. 
	start # OPTIONAL start 1 based.if not precised sample the whole VCF
	end # OPTIONAL end 1 based.if not precised sample the whole VCF
	sampling_frequency=1 # sampling frequency for site to include in the data
	missing_data="None" # Value to insert in the table in case of missing data
	bgzip=True # the file is bgzip
	write=True # write an output file? else return a table
	outfile="default" # path to the output file. if not specified create a meaningfule name :
	outfile="/".join(vcf_file.split("/")[:-1])+"/"+os.path.splitext(vcf_file.split("/")[-1])[0]+"_"+feature+"_"+chrom+"_"+str(start)+"_"+str(end)+".txt"
	

	>>> extract_features(vcf_file=sample_vcfgz,feature="DP",chrom="N00180",start=1,end=10,sampling_frequency=1,missing_data="None",bgzip=True,write=False,outfile="")[1]
	['N00180','1','2','3','None','4','10','None','None','5', '2', '1', '12', '7', 'None', '1', '5', 'None', '6', '1', '2']
	>>> extract_features(vcf_file=sample_vcfgz,feature="DP",chrom="N00180",start=1,end=10,sampling_frequency=2,missing_data="None",bgzip=True,write=False,outfile="")[1]#test sampling frequency and missing data
	['N00180', '2', '3', '3', '0', '6', '10', '3', '0', '5', '2', '1', '15', '7', '0', '1', '6', '0', '6', '2', '2']
	>>> extract_features(vcf_file=sample_vcfgz,feature="GT",chrom="N00180",start=1,end=10,sampling_frequency=10,missing_data="./.",bgzip=True,write=False,outfile="")[1]#test another feature
	['N00180', '10', '0/1', '1/1', './.', '0/1', '0/1', '0/1', '1/1', '0/1', '0/1', '0/1', '1/1', '0/1', '0/0', '1/1', '0/1', '0/1', '1/1', '0/1', '0/1']
	"""
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)
	#initialize big_table
	header=["#Chrom","Pos"]+  input_vcf.samples
	if write==False:
		feature_table=[]
		feature_table.append(header)
	#fill the big table
	elif write==True:
		if outfile=="default":
			outfile="/".join(vcf_file.split("/")[:-1])+"/"+os.path.splitext(vcf_file.split("/")[-1])[0]+"_"+feature+"_"+chrom+"_"+str(start)+"_"+str(end)+".txt"#modify default name to something meaningful in the vcf directory
		output=open(outfile,"a")
		output.write("\t".join(header)+"\n")
	else:
		raise Exception("invalid write argument")
	i=0#iterator for sampling frequency
	for record in input_vcf.fetch(chrom,start,end):# for every site
		i+=1
		if 1%100000==0: print i
		if i%sampling_frequency==0:#check if we are at sampling frequency
			feats=[str(record.CHROM.split("chr")[1]),str(record.POS)] #initialize the line with seq and position
			if feature in record.FORMAT:#check if the feature is in the line format, otherwise return missing data for every ind
				#print record.samples
				for sample in record.samples:
					if sample[feature]!= None:# if the feature exist for this ind
						feats.append(str(sample[feature]))
					else: #if this nd is missing return missing data
						feats.append(str(missing_data))			
			else : # if the feature is not in the line format return missing data for every ind
				feats=feats+list(numpy.repeat(missing_data,len(record.samples)))
			if write==True:
				output.write("\t".join(feats)+"\n")
			else:
				feature_table.append(feats)#append the line to the big tables)		
	if write==False:
		return feature_table
	print "done"


def extract_pi_region(vcf_file,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,min_nsites=0,min_variants=0,verbose="min",called=True,output="pi"):
	"""return the pi value for a region in a vcf provided a number of conditions. Or N
	One can specify which individuals he is interested in and a coverage range in which where individual should be to calculate pi.
	Finally one can define a minimum number of sites that must respect criteria to be able to compute pi as well (if not respected return None)
	as a minimum number of variants.

	vcf_file(str)the path to the vc_file of interest
	###region
	chrom(str)# chromosome of interest
	start(int)# start 1 based
	end(int)# end 1 based
	###coverage range in which every selected individual should be
	mincov(int)# inferior limit of the coverage that every individual should have
	maxcov(int)(## superior limit of the coverage that every individual should have
	called(bool)	#all individuals need an assigned genotype
	output="pi"		# if output ="extended" return a list [nvariants,nsites_considered,pi]

	###number of sites that must be considered to be able to compute pi
	min_nsites(int)# min nsites that shoudl respect criteria to ba able to return pi
	min_variants(int)# max nsites that should respect criteria to be able to return pi
	##other parameters
	inds(list)# list of individuals to consider (individual names)
	bgzip(bool)#  the vcf_file is bgzipped or not3
	verbose # True|False|"min" return some info or not. by default "min" it is a single message 

	>>> extract_pi_region(vcf_file=sample_vcfgz,chrom="N00180",start=1,end=10,mincov=0,maxcov=10000,inds="all",bgzip=True,min_nsites=0,min_variants=0)
	0.048888888888888885
	>>> extract_pi_region(vcf_file,chrom,start,end,mincov=0,maxcov=100,inds="all",bgzip=True,min_nsites=0,min_variants=2)
	'NA'
	>>> extract_pi_region(vcf_file,chrom,start,end,mincov=0,maxcov=13,inds="all",bgzip=True,min_nsites=0,min_variants=0)
	0.0
	>>> extract_pi_region(vcf_file,chrom,start,end,mincov=0,maxcov=100,inds="OC_10_F",bgzip=True,min_nsites=0,min_variants=0)# 10 sites 1 snps
	0.1
	"""
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	if inds=="all" or inds==["all"]:inds=input_vcf.samples# transform "all" in a list of all individuals in the vcf
	#Function
	pi_values=[]#list 
	nsites_considered=0#iterator for sampling frequency
	total_nsites=0
	nvariants=0# iterator for sites that are varying
	###identify individual to remove when calculating stats
	inds_to_delete=[]
	for i,ind in enumerate(input_vcf.samples):#check which ind is ion sample and compare it to our list of inds
		 if ind not in inds:#delete this ind
		 	inds_to_delete.append(i)
	#go along the region
	if chrom!="all":
		for record in input_vcf.fetch(chrom,start,end):# for every site
			cond=checkRecord_Cov(input_vcf,record,mincov,maxcov,inds=inds,called=True,nalleles=[1,2])# check if the site respect our condition
			total_nsites+=1
			if cond:# if it does
				nsites_considered+=1 
			 	if total_nsites%100000==0: print total_nsites,"sites",nsites_considered,"sites passed filter"
			 	for index in  sorted(inds_to_delete)[::-1]:#remove the individuals we do not want
			 		del record.samples[index]
			 	if verbose==True:print record.POS
			 	if verbose==True:print "inds",inds		 	
			 	if verbose==True:print "GT",[sample["GT"] for  sample in record.samples] 
			 	if verbose==True:print "DP",[sample["DP"] for  sample in record.samples]
				pi_values.append(record.nucl_diversity)#calculate pi
				if record.nucl_diversity>0.0:nvariants+=1
			#compute total information for the window
	elif chrom=="all":
		for record in input_vcf:# for every site
			cond=checkRecord_Cov(input_vcf,record,mincov,maxcov,inds=inds,called=True,nalleles=[1,2])# check if the site respect our condition
			total_nsites+=1
			if cond:# if it does
				nsites_considered+=1
			 	if total_nsites%100000==0: print total_nsites,"sites",nsites_considered,"sites passed filter"
			 	for index in  sorted(inds_to_delete)[::-1]:#remove the individuals we do not want
			 		del record.samples[index]
			 	if verbose==True:print record.POS
			 	if verbose==True:print "inds",inds		 	
			 	if verbose==True:print "GT",[sample["GT"] for  sample in record.samples] 
			 	if verbose==True:print "DP",[sample["DP"] for  sample in record.samples]
				pi_values.append(record.nucl_diversity)#calculate pi
				if record.nucl_diversity>0.0:nvariants+=1
	if verbose==True or verbose=="min":print "nvariants:",nvariants,"nsites_considered:",nsites_considered
	if output=="pi":
		if nsites_considered>=min_nsites and nvariants>=min_variants and len(pi_values):
			pi_value=sum(pi_values)/nsites_considered		
			return pi_value
		else:
			return "NA"
	elif output=="extended":
		if nsites_considered>=min_nsites and nvariants>=min_variants and len(pi_values):
			pi_value=sum(pi_values)/nsites_considered		
			return [nvariants,nsites_considered,pi_value]
		else:
			return [nvariants,nsites_considered,"NA"]
	else:
		raise Exception("incorrect output argumnent, should be pi or extended")

def extract_fasta_region(vcf_file,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,variants="N",missing_char="N"):
	"""return a dict of sequences for the region in a vcf for the individuals of interest. 
	One can specify which individuals he is interested in and a coverage range in which bases should be considered. It is possible to have a missing 
	data character that vary.
	
	For variants, one can choose to keep both variants (DIP), one of the two allele at random(RANDOM). If any other string is specified it will 
	be the character for variants.
	vcf_file(str)the path to the vc_file of interest

	vcf_file # the vcf_file from which to extract the sequence
	###region of interest
	chrom # the sequence name (scaffold or chrom)
	start # the start of the region, 1 based
	end # the end of the region, inclusive
	#conditions
	mincov=0 # the minimum coverage for a site in an individual not too be considered missing data
	maxcov=10000 # the maximum coverage for a site in an individual not too be considered missing data
	inds="all" # a list of all individuals for which you want to have the sequence in this file.If "all". All the individuals are computed
	bgzip=True #Wether the input vcf is bgzipped or not
	verbose=True
	#Variants
	variants="N" #For variants, one can choose to keep both variants ("DIP"), one of the two allele at random(RANDOM). If any other string is specified it will 
				  #be the character for variants.
	missing_char="N" # The character for a site that do not respect the coverage range or that is not called

	>>> extract_fasta_region(vcf_file=sample_vcfgz,chrom="N00180",start=1,end=11,mincov=0,maxcov=10000,inds="all",bgzip=True,variants="P",missing_char="N")
	{'>N00180_1_11_OC_8_F': 'ACAAAAAAAPC', '>N00180_1_11_OC_4_M': 'ACAAAAAAAPC', '>N00180_1_11_OC_HB9_F': 'ACAAAAAAAPP', '>N00180_1_11_OC_6_F': 'NNNNNNNAAAC', '>N00180_1_11_OC_HB2_M': 'ACAAAAAAAPP', '>N00180_1_11_OC_5_M': 'NCAAAAAAAPC', '>N00180_1_11_OC_HB7_F': 'ACAAAAAAAAC', '>N00180_1_11_OC_HB8_F': 'ACAAAAAAAPC', '>N00180_1_11_OC_2_M': 'NNNNNNNNNNN', '>N00180_1_11_OC_HB4_M': 'ACAAAAAAAAC', '>N00180_1_11_OC_HB3_M': 'NNNAAAAAACP', '>N00180_1_11_OC_HB6_F': 'NNAAAAAAAPC', '>N00180_1_11_OC_3_M': 'ACAAAAAAAPP', '>N00180_1_11_OC_HB5_M': 'ACAAAAAAAPC', '>N00180_1_11_OC_HB1_M': 'ACAAAAAAAAP', '>N00180_1_11_OC_10_F': 'ACAAAAAAAPC', '>N00180_1_11_OC_1_M': 'ACAAAAAAAAC', '>N00180_1_11_OC_7_F': 'ACAAAAAAAPC', '>N00180_1_11_OC_HB10_F': 'ACAAAAAAAPC'}
	>>> extract_fasta_region(vcf_file=sample_vcfgz,chrom="N00180",start=1,end=11,mincov=0,maxcov=10000,inds=["OC_HB6_F"],bgzip=True,variants="DIP",missing_char="N")
	{'>N00180_1_11_OC_HB6_F_allele1': 'NNAAAAAAACC', '>N00180_1_11_OC_HB6_F_allele2': 'NNAAAAAAAAC'}
	"""
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	if inds=="all" or inds==["all"]:inds=input_vcf.samples# transform "all" in a list of all individuals in the vcf
	if type(inds) == str: inds=[inds]
	dict_seq={}#dictionnary to stock diploid seq
	if variants=="DIP":
		for ind in inds:
			dict_seq[ind] = ["",""]
	else:
		for ind in inds:
			dict_seq[ind] = ""
	if not all(ind  in input_vcf.samples for ind in inds): raise Exception("not all the individuals in",inds, " are found in the vcf samples:",input_vcf.samples) 
	#Function
	###identify individual to remove when calculating stats
	inds_to_delete=[]
	for i,ind in enumerate(input_vcf.samples):#check which ind is in sample and compare it to our list of inds
		 if ind not in inds:#delete this ind
		 	inds_to_delete.append(i)
	#go along the region
	for record in input_vcf.fetch(chrom,start,end):# for every site
	 	for index in  sorted(inds_to_delete)[::-1]:#remove the individuals we do not want
	 		del record.samples[index]
	 	if "DP" in record.FORMAT:
			for sample in record.samples:
				if mincov<sample["DP"]<maxcov and sample.called==True:
					if variants=="DIP":
						dict_seq[sample.sample][0]+=sample.gt_bases.split("/")[0]
						dict_seq[sample.sample][1]+=sample.gt_bases.split("/")[1]
					elif variants=="RAN":#randomly pick allele one or two every time
						dict_seq[sample.sample]+=sample.gt_bases.split("/")[random.choice([0,1])]
					else :
						if sample.gt_bases.split("/")[0]!=sample.gt_bases.split("/")[1]: # If the two alleles are different add the character specify in "variants" 
							dict_seq[sample.sample]+=variants
						else:
							dict_seq[sample.sample]+=sample.gt_bases.split("/")[0]
				else:
					if variants=="DIP":
						dict_seq[sample.sample][0]+=missing_char
						dict_seq[sample.sample][1]+=missing_char
					else:
							dict_seq[sample.sample]+=missing_char
		else:
			if variants=="DIP":
				for key in dict_seq.keys():
					dict_seq[key][0]+=missing_char
					dict_seq[key][1]+=missing_char
			else:
				for key in dict_seq.keys():
					dict_seq[key]+=missing_char
	#Cheange the key to fasta header
	final_dict={}
	for key in dict_seq.keys():
		newkey=">"+chrom+"_"+str(start)+"_"+str(end)+"_"+key
		if variants=="DIP":
			newkey1=">"+chrom+"_"+str(start)+"_"+str(end)+"_"+key+"_allele1"
			final_dict[newkey1]=dict_seq[key][0]
			newkey2=">"+chrom+"_"+str(start)+"_"+str(end)+"_"+key+"_allele2"
			final_dict[newkey2]=dict_seq[key][1]
		else:
			newkey=">"+chrom+"_"+str(start)+"_"+str(end)+"_"+key
			final_dict[newkey]=dict_seq[key]
	return final_dict

def sum_pairwise_differences(vcf_file,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,called=True,output="sum",nb_ind_with_min_cov="all"):
	"""return the pi value for a region in a vcf provided a number of conditions. 
	One can specify which individuals he is interested in and a coverage range in which where individual should be to calculate pi.
	Finally one can define a minimum number of sites that must respect criteria to be able to compute pi as well (if not respected return None)
	as a minimum number of variants.

	vcf_file(str)the path to the vc_file of interest
	###regiom
	chrom(str)# chromosome of interest
	start(int)# start 1 based
	end(int)# end 1 based
	###coverage range in which every selected individual should be
	mincov(int)# inferior limit of the coverage that every individual should have
	maxcov(int)(## superior limit of the coverage that every individual should have
	called(bool)	#all individuals need an assigned genotype
	output="sum"		# if output ="extended" return a list [nsitesok,,pi]

	###number of sites that must be considered to be able to compute pi
	min_nsites(int)# min nsites that shoudl respect criteria to ba able to return pi
	min_variants(int)# max nsites that should respect criteria to be able to return pi
	##other parameters
	inds(list)# list of individuals to consider (individual names)
	bgzip(bool)#  the vcf_file is bgzipped or not3
	verbose # True|False|"min" return some info or not. by default "min" it is a single message 
	"""
	###CHOOSE THE RIGHT VCF
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	if inds=="all" or inds==["all"]:inds=input_vcf.samples# transform "all" in a list of all individuals in the vcf
	#Function
	sum_pairwise=0#iterator for sampling frequency
	nsites_ok=0
	###identify individual to remove when calculating stats
	inds_to_delete=[]
	for i,ind in enumerate(input_vcf.samples):#check which ind is ion sample and compare it to our list of inds
		 if ind not in inds:#delete this ind
		 	inds_to_delete.append(i)
	#go along the region
	if chrom!="all":
		check=len(sh.tabix(vcf_file,str(chrom)+":"+str(start)+"-"+str(end)))
		#print "check;' ",check,"'"
		if check==0: 
			if output=="sum":
				return 0
			elif output=="extended":
				return [0,0]
		for record in input_vcf.fetch(chrom,start,end):# for every site
			#print "HERE"
			cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=[1,2],nb_ind_with_min_cov=nb_ind_with_min_cov)# check if the site respect our condition
			#print "cond",cond
			#print "HERE2"
			if cond:# if it does
			 	for index in  sorted(inds_to_delete)[::-1]:#remove the individuals we do not want
			 		del record.samples[index]
			 	#print record.samples
				nsites_ok+=1
				#print record.nucl_diversity
				sum_pairwise+=record.nucl_diversity 
			#compute total information for the window
	elif chrom=="all":
		for record in input_vcf:# for every site
			cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=[1,2],nb_ind_with_min_cov=nb_ind_with_min_cov)# check if the site respect our condition
			if cond:# if it does
			 	for index in  sorted(inds_to_delete)[::-1]:#remove the individuals we do not want
			 		del record.samples[index]
			 	#print record.samples
				nsites_ok+=1
				sum_pairwise+=record.nucl_diversity 
			#compute total information for the window
	if output=="sum":
		return sum_pairwise
	elif output=="extended":
		return [sum_pairwise,nsites_ok]
	#Go in normal vcf and count sites
	
def count_sites_under_condition_vcf(vcf_file,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,nb_ind_with_min_cov="all",nalleles=[1,2],snps=False):
	"""
	count the number of sites in a vcf file in a given stretch that respect a certain condition
	output a list [nsites_pass_the_filter, nsites_total_in_the_stretch_considered]

	vcf_file #vcf_file
	chrom="all"|scaf (str)# can be all and then end and start are disregarded the whole vcf is considered  
	start# the start position (1 based)
	end# the end position, inclusive
	mincov=0# the min coverage for every ind to consider a site okay
	maxcov=10000# he max coverage for every ind to consider a site okay
	inds # a list of selected individuals. if inds="all" or ["alls"] consider every sample
	bgzip=True# is the file bgzipped
	nb_ind_with_min_cov="all" # the number of individuals that need at least mincov to call the site. If "all", it becomes the same as ind, If "all_vcf", it becomes all the individuals in the vcf  Not applicable to maxcov!
	nalleles = [1,2]  the number of alleles the sites can have, in this case one or two, it is used to determine the amount of true variable sites out of snps positions
	>>> count_sites_under_condition_vcf(sample_vcfgz,chrom="N00180",start=100,end=199,mincov=0,maxcov=10000,inds="all",bgzip=True)
	N00180 100 199
	[100, 100]
	>>> count_sites_under_condition_vcf(sample_vcfgz,chrom="N00180",start=100,end=199,mincov=7,maxcov=10000,inds="all",bgzip=True)
	N00180 100 199
	[29, 100]
	snps , if true, only consider sites that are also snps
	"""
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	nsites_OK=0
	nsites_total=0
	#print "in count_sites_under_condition_vcf nb_ind_with_min_cov :",nb_ind_with_min_cov, " inds", ind
	if chrom!="all":
			print chrom,start,end
			check=len(sh.tabix(vcf_file,str(chrom)+":"+str(start)+"-"+str(end)))
			#print  check
			#print "check;' ",check,"'"
			if check==0: 
				return [0,0]
			for record in input_vcf.fetch(chrom,start,end):# for every site
				cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=nalleles,nb_ind_with_min_cov=nb_ind_with_min_cov,snps=snps)# check if the site respect our condition
				nsites_total+=1
				if cond:# if it does
					#if  any([int(sample['DP'])<5 for sample in record.samples]): print [int(sample['DP']) for sample in record.samples] # to check this argument nb_ind_with_min_cov
					nsites_OK+=1
	elif chrom=="all":
		for record in input_vcf:# for every site
			cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=nalleles,nb_ind_with_min_cov=nb_ind_with_min_cov,snps=snps)# check if the site respect our condition
			nsites_total+=1
			if cond:# if it does
				nsites_OK+=1
	return [nsites_OK,nsites_total]

def count_sites_under_condition_vcf_allbed(input_bed,input_vcf,output_file,mincov=0,maxcov=100000000,inds="all",nalleles=[1,2,3,4],bgzip=True):
	'''	count the number of sites in a vcf file across a whole bed that respect a certain condition
	output a list [nsites_pass_the_filter, nsites_total_in_the_stretch_considered]

	input_bed # the bed file that contains all the regions over which you wanna sum the sites that pass a given filter
	input_vcf # the vcf with the sited
	output_file # the output file
	mincov=0# the min coverage for every ind to consider a site okay
	maxcov=10000# he max coverage for every ind to consider a site okay
	inds # a list of selected individuals. if inds="all" or ["alls"] consider every sample
	nalleles = [1,2,3,4]  # a list with the number of alleles a site can have to be considered for example nalleles = [1,2] means that sites with 3 or 4 alleles are discarded
	bgzip=True# is the file bgzipped
	'''
	output=open(output_file,"w")
	with open(input_bed) as f :
		for line in f:
			line_split=line.split()
			nsites_ok , nsites_total = count_sites_under_condition_vcf(input_vcf,line_split[0],int(line_split[1])+1,int(line_split[2]),mincov=mincov,maxcov=maxcov,inds="all",bgzip=True)
			line_split.append(str(nsites_ok))
			line_split.append( str(nsites_total))
			print line_split
			output.write("\t".join(line_split)+"\n")
	output.close()

def checkSnp_Cov(input_vcf,record,mincov=0,maxcov=100000000,inds="all",nalleles=[1,2,3,4],nb_ind_with_min_cov="all",snps=False):
	"""Check if a site respect coverage conditions for selected individuals. Return TRUE/False. In order to be true if any coverage is specified the individual has to be called

		input_vcf # an object of class vcf.parser.Reader from which the recorded is extracted
		record #  an object of class  vcf.model._Record
		mincov # the coverage minimum the site should have for the for every  individual sppecified with ind
		maxcov # the coverage maximum the site should have for the for every  individual sppecified with ind
		inds # a list of selected individuals. if inds="all" or ["alls"] consider every sample
		called	#all individuals need an assigned genotype
		nalleles (list) # a list of integer of alleles tolerated for the check. for example [1] is for monomorphic snp, [1,2] include snps and [1,2,3] include triallele snps
		nb_ind_with_min_cov="all" # the number of individuals that need at least mincov to call the site. If "all", it becomes the same as ind, If "all_vcf", it becomes all the individuals in the vcf  Not applicable to maxcov!
		>>> input_vcf=vcf.Reader(fsock=None, filename=sample_vcfgz, compressed=True, prepend_chr="False", strict_whitespace=False)
		>>> record=input_vcf.next()
		>>> checkSnp_Cov(input_vcf,record)###check function
		True
		>>> checkSnp_Cov(input_vcf,record,5)#check mincov
		False
		>>> checkSnp_Cov(input_vcf,record,maxcov=12)# checkmaxcov
		True
		>>> checkSnp_Cov(input_vcf,record,maxcov=5)#check maxcov
		False
		>>> checkSnp_Cov(input_vcf,record,maxcov=10,inds=['OC_3_M','OC_4_M','OC_5_M'])#check inds
		True
		>>> checkSnp_Cov(input_vcf,record,5,nb_ind_with_min_cov=3)
		True
	"""
	###Checks
	#print "in checkSnp_Cov nb_ind_with_min_cov :",nb_ind_with_min_cov, " inds", inds

	if type(input_vcf)!=vcf.parser.Reader: raise Exception ("input_vcf must be a parser.Reader object")
	if type(record)!=vcf.model._Record: raise Exception ("record must be a  vcf.model._Record object")
	#function
	if inds=="all" or inds==["all"]:inds=input_vcf.samples# if no list of individuals, us all individuals 
	if nb_ind_with_min_cov=="all" : nb_ind_with_min_cov=len(inds)# if we want all the individuals with at least X coverage
	if nb_ind_with_min_cov=="all_vcf" : 
		nb_ind_with_min_cov=len(input_vcf.samples)# if we want all the individuals with at least X coverage
		inds=input_vcf.samples
	#print "in checkSnp_Cov nb_ind_with_min_cov :",nb_ind_with_min_cov, " inds", inds
	if not len(record.alleles) in nalleles: return False # check the number of alleles
	#print snps 
	if snps==True and record.num_het == 0: return False
	if "DP" in record.FORMAT:
		if mincov==0 :# we want to avoid to take in the None that are 0 and prevent us to use the condition
			cond=all([ ind["DP"]<=maxcov for ind in record.samples if ind["DP"]!=None and (ind.sample in inds) ]) 
		else: # we need to take into account the none cause they mean DP=0 and cond=False
			#print [mincov<= ind["DP"]<=maxcov and ind.called==True for ind in record.samples  if (ind.sample in inds)].count(True)   
			cond=nb_ind_with_min_cov<=[mincov<= ind["DP"]<=maxcov and ind.called==True for ind in record.samples  if (ind.sample in inds)].count(True)   
	else:
		cond=False#raise Exception ("format do not include DP for",record.POS,record.CHROM)
	#if cond==True: print [sample["DP"] for sample in record.samples if (sample.sample in inds)]
	#if cond==True: assert 8<=[5<= ind["DP"]<=maxcov and ind.called==True for ind in record.samples  if (ind.sample in inds)].count(True)," very ugly.... just to check that 1 individual is okay for 7 reads and that I check that all along"
	return cond

def pi_double_vcf(vcf_file_snps,vcf_allsites,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,min_nsites=0,max_nsites=10000000000,called=True,nb_ind_with_min_cov="all"):
	''' calculate pi of a region using two vcfs, one for SNPs and one for allsites ( could be the same vcf)

	vcf_file_snps # the path the the snps files
	vcf_allsites # the path to the allsites files
	chrom # sequence ( chromosome /scaffold) of the region
	start # start of the region
	end # end
	mincov=0 # 
	maxcov=10000
	inds="all"
	bgzip=True
	min_nsites=0 # the minimum number of valid sites in any given window for pi to be estimated, if less return NA
	max_nsites=10000000000  # the maximum number of valid sites in any given window for pi to be estimated, if less return NA
	called=True# all individuals should be assigned a genotype (for SNPs)
	nb_ind_with_min_cov="all"  # the number of individuals in "inds" that should have the mincov for a given site, by default all individuals#

	'''
	#print "pi_doublevcf() without subsampling"
	count_sites = count_sites_under_condition_vcf(vcf_allsites,chrom,start,end,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov)
	count_snps = count_sites_under_condition_vcf(vcf_file_snps,chrom,start,end,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov,nalleles=[2],snps=True) # this line just serve as giving an extra info in the logs files
	#print "pi_double_vcf()sum pairwise"
	if count_sites[0] >max_nsites:
		a = subsample_to_X_callable_sites(vcf_allsites,chrom,start,end,max_nsites,mincov,maxcov,inds)
		sum_pairwise = sum_pairwise_differences(vcf_file_snps,a[0],a[1],a[2],mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov)
		count_sites = count_sites_under_condition_vcf(vcf_allsites,a[0],a[1],a[2],mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov)
		count_snps = count_sites_under_condition_vcf(vcf_file_snps,a[0],a[1],a[2],mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov,nalleles=[2],snps=True) # this line just serve as giving an extra info in the logs files

	else:
		#print "enter 1"
		sum_pairwise = sum_pairwise_differences(vcf_file_snps,chrom,start,end,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov)
	if max_nsites >= count_sites[0]>=min_nsites  :
		#print "enter 2"
		pi = sum_pairwise/count_sites[0]
	else:
		pi = "NA"
	print "pi_double_vcf()",pi,"nsites_ok",count_sites[0],"win_len",count_sites[1] ,end-start
	return [pi,count_sites[0],count_sites[1],count_snps[0]] # pi, nsites passing condition, nsites, nsnps

def subsample_to_X_callable_sites(vcf_file,chrom,start,end,max_nsites,mincov,maxcov,inds,bgzip=True,nb_ind_with_min_cov="all"):
	''' subsample a region to the smallest region from the start that can satisfy the conditions requested and output the redyced region ( maximum number of callable sites)
	
	vcf_file = the vcf file
	chrom # the sequence (chromosome /scaffold) where the unsampled region is
	start # the start of the unsampled region
	end # the end of the unsampled region
	max_nsites # the maximum number of valid sites for the regions
	mincov #the minimum coverage for a given site that all individuals must have for the site to be considered but see  nb_ind_with_min_cov )
	maxcov#the maximum coverage, if any individual in "inds" have more, the site is excluded
	inds # a list of individuals you want to estimate nucleotide diversity for,
	bgzip  #wether the vcffiles are gzipped or not
	nb_ind_with_min_cov="all"
	'''
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	nsites_OK=0
	nsites_total=0
	first_site=[]
	for record in input_vcf.fetch(chrom,start+1,end):# for every site
		cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=[1,2],nb_ind_with_min_cov=nb_ind_with_min_cov)# check if the site respect our condition
		nsites_total+=1
		if cond:# if it does
			nsites_OK+=1
			if not first_site:
				first_site=record.POS
		if nsites_OK == max_nsites:
			print "subsample()",chrom,first_site,record.POS
			return [chrom,first_site,record.POS]
	else:
		print nsites_OK
		return None # no possible subsampling


def calculate_pi_of_a_bed_double_vcf(bed,vcfsnps,vcf_allsites,mincov,maxcov,bgzip,inds,output_file):
	''' calculate the average pi over all windows in a bed 
	
	bed # the bed containing the window for which you want the sites
	vcfsnps # the path the the snps files
	vcf_allsites # the path to the allsites files
	mincov #the minimum coverage for a given site that all individuals must have for the site to be considered but see  nb_ind_with_min_cov )
	maxcov#the maximum coverage, if any individual in "inds" have more, the site is excluded
	bgzip  #wether the vcffiles are gzipped or not
	inds # a list of individuals you want to estimate nucleotide diversity for,
	output_file # the output file
	'''
	nsitesok = 0
	sum_pairwise = 0
	with open (bed) as f:
		for line in f:
			print line
			sum_pairwise+= sum_pairwise_differences(vcfsnps,chrom=line.split()[0],start=int(line.split()[1])+1,end=int(line.split()[2]),mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip)
			nsitesok += count_sites_under_condition_vcf(vcf_allsites,chrom=line.split()[0],start=int(line.split()[1])+1,end=int(line.split()[2]),mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip)[0]
			print sum_pairwise,nsitesok
	if nsitesok==0: 
		pi=0
	else:
		pi= sum_pairwise/nsitesok
	output=open(output_file,"a")
	output.write(bed+"\t"+str(pi)+"\t"+str(nsitesok)+"\n")
	output.close()



def pi_double_vcf_fix_nsites(vcf_file_snps,vcf_allsites,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,min_nsites=0,max_nsites=10000000000,called=True):
	'''
	'''
	#print "pi_doublevcf() without subsampling"
	#count_sites = count_sites_under_condition_vcf(vcf_allsites,chrom,start,end,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip)
	#print "pi_double_vcf()sum pairwise"
	print vcf_allsites,chrom,start,end,max_nsites,mincov,maxcov,inds
	a = subsample_to_X_callable_sites(vcf_allsites,chrom,start,end,max_nsites,mincov,maxcov,inds)
	print a 
	sum_pairwise = sum_pairwise_differences(vcf_file_snps,a[0],a[1],a[2],mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip)
	print a
	pi = sum_pairwise/int(max_nsites)
	print "pi_double_vcf()",pi,max_nsites,"NA"
	return [pi,max_nsites,"NA"] # pi, nsites passing condition, total nsites


def filtervcf(vcf_file,vcf_output,mincov=1,maxcov=10000,bgzip =True,inds="all",nb_ind_with_min_cov="all",nalleles=[1,2],snps=False):
	''' filter a gzipped vcf'''
	output_handle = gzip.open(vcf_output,"wb")
	f = gzip.open(vcf_file)
	for line in f:
		if line.startswith("#"):
			output_handle.write(line)
		else:
			break
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	for record in input_vcf:# for every site
		cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=nalleles,nb_ind_with_min_cov=nb_ind_with_min_cov,snps=snps)# check if the site respect our condition
		if cond:
			output_handle.write(line)
			line = f.readline()
		else:
			line = f.readline()
	output_handle.close()


def create_vcf_from_annot(vcfinput,bed,annot,vcfoutput):
	'''create a subset of vcfinput that contains only the sites mentionned
	 in the bedfile associated to a given annotation in that file'''
	# get the header
	output_handle = gzip.open(vcfoutput,"wb")
	with gzip.open(vcfinput) as f:
		for line in f:
			if line.startswith("#"):
				output_handle.write(line)
			else:
				break
	output_handle.close()
	with open(bed) as f:
		for line in f:
			if annot in line:
				os.system("tabix "+vcfinput+" "+line.split()[0]+":"+str(int(line.split()[1])+1)+"-"+str(line.split()[2])+" | gzip -c >> "+vcfoutput) 




def count_sites_under_condition_vcf_to_set(vcf_file,chrom,start,end,mincov=0,maxcov=10000,inds="all",bgzip=True,nb_ind_with_min_cov="all",nalleles=[1,2],snps=False):
	"""
	return a set of a ll the sites in a vcf file that respect a certain condition that respect a certain condition
	output a set

	vcf_file #vcf_file
	chrom="all"|scaf (str)# can be all and then end and start are disregarded the whole vcf is considered  
	start# the start position (1 based)
	end# the end position, inclusive
	mincov=0# the min coverage for every ind to consider a site okay
	maxcov=10000# he max coverage for every ind to consider a site okay
	inds # a list of selected individuals. if inds="all" or ["alls"] consider every sample
	bgzip=True# is the file bgzipped
	nb_ind_with_min_cov="all" # the number of individuals that need at least mincov to call the site. If "all", it becomes the same as ind, If "all_vcf", it becomes all the individuals in the vcf  Not applicable to maxcov!
	nalleles = [1,2]  the number of alleles the sites can have, in this case one or two, it is used to determine the amount of true variable sites out of snps positions
	"""
	set_ok_sites = set()
	input_vcf=vcf.Reader(fsock=None, filename=vcf_file, compressed=bgzip, prepend_chr="False", strict_whitespace=False)#open the vcf parser
	nsites_OK=0
	nsites_total=0
	#print "in count_sites_under_condition_vcf nb_ind_with_min_cov :",nb_ind_with_min_cov, " inds", ind
	if chrom!="all":
			#print chrom,start,end
			check=len(sh.tabix(vcf_file,str(chrom)+":"+str(start)+"-"+str(end)))
			#print  check
			#print "check;' ",check,"'"
			if check==0: 
				return [0,0]
			for record in input_vcf.fetch(chrom,start,end):# for every site
				cond=checkSnp_Cov(input_vcf,record,mincov,maxcov,inds=inds,nalleles=nalleles,nb_ind_with_min_cov=nb_ind_with_min_cov,snps=snps)# check if the site respect our condition
				nsites_total+=1
				if cond:# if it does
					#if  any([int(sample['DP'])<5 for sample in record.samples]): print [int(sample['DP']) for sample in record.samples] # to check this argument nb_ind_with_min_cov
					set_ok_sites.add(str(record.CHROM)+"_"+str(record.POS))
	return set_ok_sites

def make_subset_bed_annot(input_bed,annot,output_bed):
	''' create an output_bed file that is the subset of the input_bed
	file that contains only the line that match the pattern  annot
	'''
	output = open(output_bed,"w")
	with open(input_bed) as f:
		for line in f:
			if annot in line:
				output.write(line)
	output.close()

def sum_pairwise_over_bed(vcf_file,bed,mincov=1,maxcov=10000,inds="all",bgzip=True,called=True,nb_ind_with_min_cov="all"):
	''' function sum_pairwise_differences( see vcf_tools documentation) iterated over all regions in a bed file
	return the sum of pairwise differences, and the number of sites considered
	'''
	sum_overall,nsites = 0,0
	with open(bed) as f:
		for line in f:
			print int(line.split()[2])-int(line.split()[1]),line
			seq,start,end = line.split()[0],int(int(line.split()[1])+1),int(line.split()[2])
			temp = sum_pairwise_differences(vcf_file,chrom=seq,start = start,end = end ,mincov=mincov,maxcov=maxcov,inds= inds,bgzip=bgzip,called=called,output="extended",nb_ind_with_min_cov=nb_ind_with_min_cov)
			sum_overall += temp [0]
			nsites += temp[1]
	return sum_overall,nsites




def sites_under_condition_vcf_from_bed_to_set(bed,vcf_file,mincov=0,maxcov=10000,inds="all",bgzip=True,nb_ind_with_min_cov="all",nalleles=[1,2],snps=False):
	''' a version of count_sites_under_condition_vcf that return a set for all sites respecting a condition that intersects with a bed'''
	set_sites_ok = set()
	with open(bed) as f:
		for line in f:
			seq,start,end = line.split()[0],int(int(line.split()[1])+1),int(line.split()[2])
			temp_set =  count_sites_under_condition_vcf_to_set(vcf_file,seq,start,end,mincov=mincov,maxcov=maxcov,inds = inds,bgzip=bgzip,nb_ind_with_min_cov= nb_ind_with_min_cov ,nalleles= nalleles ,snps = snps)
			new_set = set_sites_ok.union(temp_set) # add all the other sites
			set_sites_ok = new_set
	return set_sites_ok



def condition_overlap_over_Multiplevcfs(bed,vcf_list =[],mincov=1,maxcov=10000,inds="all",bgzip=True,called=True,nb_ind_with_min_cov="all",nalleles= [1,2],snps=False):
	''' a version of count_sites_under_condition_vcf that return a set for all sites respecting a condition across several vcfs  for regions defined in a  bed.
	'''
	set_sites_ok = sites_under_condition_vcf_from_bed_to_set(bed,vcf_list[0],mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov,nalleles=nalleles,snps=snps)
	print vcf_list[0], len(set_sites_ok)
	for vcf_file in vcf_list[1:]:
		temp_set = sites_under_condition_vcf_from_bed_to_set(bed,vcf_file,mincov=mincov,maxcov=maxcov,inds=inds,bgzip=bgzip,nb_ind_with_min_cov=nb_ind_with_min_cov,nalleles=nalleles,snps=snps) 
		set_sites_ok = set_sites_ok.intersection(temp_set)
		print vcf_file, len(temp_set)
	return set_sites_ok

#Special function used for SBE project


def extract_pi_double_vcf_bed_for_several_beds_and_several_allsites_files(list_beds,vcf_snps,vcf_allsites_list,output_file,annot="",mincov=3,maxcov=10000,bgzip=True,called=True,nb_ind_with_min_cov="all"):
	'''A Special function that is integrating many levels.  It derives from extract_pi_double_vcf_bed(). but it uses several allsites files not only 1 if the individals in the snps file are in different allsites vcfs.
	It also use a lit file  list_beds that contains a lit of bedfiles that one want to caculate vcf for. In each of this file, only the line containing the string defined bu the parameter annot will be used.
	NOTE: It has originally been used to calculate exonic diversity, each file in the list of file being a  gene'''
	output_handle = open(output_file,"w")
	with open(list_beds) as f:
		for filename in f:
			print filename+"\n"
			filepath = filename.strip()
			tempbed =  annot+"_"+os.path.basename(filepath)
			make_subset_bed_annot(filepath,annot,tempbed)
			print "summing"
			sums = sum_pairwise_over_bed(vcf_snps,tempbed,mincov=mincov,maxcov=maxcov,inds="all",bgzip=True,called=True,nb_ind_with_min_cov=nb_ind_with_min_cov)
			print "iterating overallsites"
			sites = condition_overlap_over_Multiplevcfs(tempbed,vcf_allsites_list ,mincov=mincov,maxcov=maxcov,inds="all",bgzip=True,called=True,nb_ind_with_min_cov=nb_ind_with_min_cov)
			#print sites
			if len(sites)>0:
				pi = sums[0]/len(sites)
			else:
				pi=  "na"
			 	sites = ""
			output_handle.write("\t".join([filepath,str(pi),str(len(sites)),"\n"]))
			os.system("rm "+tempbed)
	output_handle.close()


def test():
	pass

test()