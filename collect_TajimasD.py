#Collect TajimasD
import os
tajd_folder = "TajimasD/allgenes/"
output_table = "summary/alldata_collared.txt"

dict_taj = [] 
for filename in os.listdir(tajd_folder):
	if filename.endswith(".Tajima.D"):
		raise Exception
		gene = filename.split(".")[0]
		data = [line for line in open(tajd_folder+"/"+ filename)]
		if len(data) ==2:
			nsites =  data[1].split()[2]
			tajD =  data[1].split()[3]
		else:
			nsites = "NA"
			tajD = "NA"
		dict_taj[gene] = [nsites,tajD]

output = "temp"
with open("output_table") as f:
	for line in f:
		if line.split()[0].startswith("ENS")
			gene = ...
			newline = line.strip()+"\t"+"\t".join(dict_taj[gene])+"\n"
			output.write(newline)
		elif line.startswith("gene"):
			output.write(line.strip()+"\t"+"\tnsitestaj\ttajd".join(dict_taj[gene])+"\n")
		else:
			raise Exception

output.close()

#os.system("mv temp "+output_table)


collect the bins 