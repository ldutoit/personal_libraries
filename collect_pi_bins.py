#!/usr/bin/env python 2

import os
import vcf_tools as vcf_t
import numpy as np
oland_dict = {}
gotland_dict = {}
i=0
for filename in os.listdir("beds/CDSbeds/"):
	i+=1
	if i%100==0: print i
	if  filename.startswith("EN"):
		#gotland
		gotland_file = "pi_bigpops/gotland86allclean/"+filename.strip(".bed")+"_withpi.bed"
		if os.path.exists(gotland_file):
			gotland_data = vcf_t.get_average_pi_from_bed_output_vcf_tools(gotland_file)
			gotland_data = [gotland_data[0]*gotland_data[1],gotland_data[1]]
			gotland_dict[filename.strip(".bed")] = gotland_data
		else:
			gotland_dict[filename.strip(".bed")] = [0,0]
			print "gotland",filename
		#oland
		oland_file = "pi_bigpops/oland86clean/"+filename.strip(".bed")+"_withpi.bed"
		if os.path.exists(oland_file):
			oland_data = vcf_t.get_average_pi_from_bed_output_vcf_tools(oland_file)
			oland_data = [oland_data[0]*oland_data[1],oland_data[1]]
			oland_dict[filename.strip(".bed")] = oland_data
		else:
			oland_dict[filename.strip(".bed")] = [0,0]
			print "oland",filename#


import pickle
pickle.dump(oland_dict,open("tempdictoland.py","wb"))
pickle.dump(gotland_dict,open("tempdictgotland.py","wb"))

import pickle 
oland_dict = pickle.load(open("tempdictoland.py","rb"))
gotland_dict = pickle.load(open("tempdictgotland.py","rb"))



folder_with_gene_lists = "beds/list_genes_per_bin/"


def extract_pi_from_list_of_gene(file):
			gotland_pi =(0,0)
			oland_pi = (0,0)			
			with open(list_gene) as f:
				for gine in f:
					gene = gine.strip()
					gotland_pi += np.array(gotland_dict[gene]) 
					oland_pi += np.array(oland_dict[gene]) 
			pi = ((gotland_pi[0]/gotland_pi[1]) + (oland_pi[0]/oland_pi[1])) / 2
			return pi

#allbins
i=0
output=open("temppi1","w")
with open("summary/all_bins.txt") as f:
	for line in f:
		i+=1
		if line.startswith("sampling\t"):
			output.write(line.strip()+"\tpi\n")
		else:
			print i
			bin = line.split()[0].split("_")
			bin_number = line.split()[1]
			#extract sampling type
			sampling_type=""
			if "both" in bin: 
				sampling_type="both_absolute"
			elif "female" in bin: 
				sampling_type="female"
			else: 
				sampling_type="male"
			list_gene = '%snbins10bin%s_%s_rna_coll_%s_autosomes' % (folder_with_gene_lists,bin_number,sampling_type,bin[-2])
			print i,list_gene
			pi = extract_pi_from_list_of_gene(list_gene)
			newline = line.strip()+"\t"+str(pi)+"\n"
			output.write(newline)
output.close()


#CAREFUL 
#!mv temppi1 summary/all_bins.txt


nrep = 500
count = 0
#allbins
output=open("temppi2","w")
with open("summary/gonad_bins.txt") as f:
	for line in f:
		count+=1
		if line.startswith("sampling\t"):
			output.write(line.strip()+"\tpi\tlowCI\thighCI\n")
		else:
			bin = line.split()[0].split("_")
			bin_number = line.split()[1]
			#extract sampling type
			sampling_type=""
			if "both" in bin: 
				sampling_type="both_absolute"
			elif "female" in bin: 
				sampling_type="female"
			else: 
				sampling_type="male"
			gotland_pi =(0,0)
			oland_pi = (0,0)
			list_gene = '%snbins10bin%s_%s_rna_coll_%s_autosomes' % (folder_with_gene_lists,bin_number,sampling_type,bin[-2])
			print count,list_gene
			nonboot_pi = extract_pi_from_list_of_gene(list_gene)
			#bootstraps...
			allpis = []
			for i in range(1,nrep+1):
				list_gene =  '%snbins10bin%s_%s_rna_coll_%s_autosomesrep%i' % (folder_with_gene_lists,bin_number,sampling_type,bin[-2],i)
				pi = extract_pi_from_list_of_gene(list_gene)
				allpis.append(pi)
			lowCI,highCI = np.percentile([float(value) for value in allpis],[2.5,97.5])
			newline = line.strip()+"\t"+str(nonboot_pi)+"\t"+str(lowCI)+"\t"+str(highCI)+"\n" 
			output.write(newline)
output.close()


#CAREFUL
#!mv temppi2 summary/gonad_bins.txt
