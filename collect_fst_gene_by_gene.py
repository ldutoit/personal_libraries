import os
import vcf_tools as vcft

output_file = "summary/fst_for_anova.txt"


#get the fst gene by gene

folder_output = "fst_gene_by_gene/logs/"
full_dict = {}
for filename in os.listdir(folder_output):
	print filename
	tempdict = vcft.parce_vcftools_fst(folder_output+"/"+filename,exception=True,remove_ext = True)
	full_dict.update(tempdict) 
# just make a new dictionnary where keys are actually gene name and nom vcf files



#according to the threw list make 1 file, 2 columns --> R anova!!
output = open(output_file,"w")
output.write("gene\tmeanfst\tweightedfst\tbin_number\n")
folder_with_list = "summary/binsfst/"
for filename in os.listdir(folder_with_list):
	if "nbins3bin" in filename:
		print filename
		bin_number = filename.split("_")[0].split("bin")[2]
		with open(folder_with_list+filename) as f:
			for line in f:
				gene =line.strip()
				if  gene in full_dict.keys():
					fsts = "\t".join(full_dict[gene])
				else:
					fsts  ="NA\tNA"
				output.write(gene+"\t"+fsts+"\t"+bin_number+"\n")

output.close()

 