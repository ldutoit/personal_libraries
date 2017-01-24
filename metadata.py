#!/usr/bin/env python2
# Filename:metadata_SBE.py


##This library use different metadata I made used of throughout my PhD working with various flycatcher datasets

class SBEpop(object):
	''' a class that has file location and sex information for the project about Sex biased gene expression. 
	males9 and females 9 are the 9 individuals with the highest coverage, they are used for FST calculations while all
	individuals are used for pi
	self.ind18 is the sum of males9 and feamles 9.  '''
	def __init__(self,popcode,vcf_snps,vcf_allsites,vcf_snps18,males9,females9,allmales,allfemales,species):
		self.popcode = popcode
		self.vcf_snps = vcf_snps
		self.vcf_allsites = vcf_allsites
		self.vcf_snps18 = vcf_snps18
		self.males9 = males9
		self.females9 = females9
		self.allmales = allmales
		self.allfemales = allfemales
		self.species = species
		self.inds = allmales+allfemales
		self.inds18 = males9+females9

CZC =  SBEpop(
	popcode = "CZC",
	vcf_snps =  "/home/ludovic/nobackup_ludo/SBE/source_files/gt_CZC.rm.vcf.bgz",
	vcf_allsites = "/home/ludovic/nobackup_ludo/SBE/source_files/CZC.gatk.allsites.rm.vcf.gz",
	vcf_snps18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/gt_CZC_forfst.rm.vcf.bgz",
	males9 = ["CZC_304236_M", "CZC_382761_M", "CZC_419135_M", "CZC_83359_M", "CZC_48837_M", "CZC_82259_M", "CZC_82355_M", "CZC_82356_M", "CZC_83353_M" ],\
	females9 = ["CZC_119601_F", "CZC_382648_F", "CZC_419251_F", "CZC_91360_F", "CZC_73077_F", "CZC_82317_F", "CZC_82348_F", "CZC_82352_F", "CZC_91319_F" ],\
	allmales = ["CZC_304236_M", "CZC_382761_M", "CZC_419135_M", "CZC_83359_M", "CZC_48837_M", "CZC_82259_M", "CZC_82355_M", "CZC_82356_M", "CZC_83353_M" ,"CZC_461398_M"],\
	allfemales =["CZC_119601_F", "CZC_382648_F", "CZC_419251_F", "CZC_91360_F", "CZC_73077_F", "CZC_82317_F", "CZC_82348_F", "CZC_82352_F", "CZC_91319_F" , "CZC_419493_F"],\
	species ="collared")


I =  SBEpop(
	popcode = "I",
	vcf_snps =  "/home/ludovic/nobackup_ludo/SBE/source_files/gt_I.rm.vcf.bgz",
	vcf_allsites = "/home/ludovic/nobackup_ludo/SBE/source_files/I.gatk.allsites.rm.vcf.gz",
	vcf_snps18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/gt_I_forfst.rm.vcf.bgz",
	males9 = ["I_D5-2_M", "I_D6-2_M", "I_DA-3_M", "I_DB-1_M", "I_E30-3_M", "I_LC-2_M", "I_LX-2_M", "I_LZ-1_M", "I_S1_M" ],\
	females9 = ["I_F1-1_F", "I_F16-2_F", "I_F9-1_F", "I_G13-2_F", "I_GX_F", "I_S4-19_M", "I_LA_F", "I_LY-1_F", "I_S3_F" ],\
	allmales = ["I_D5-2_M", "I_D6-2_M", "I_DA-3_M", "I_DB-1_M", "I_E30-3_M", "I_LC-2_M", "I_LX-2_M", "I_LZ-1_M", "I_S1_M" ],\
	allfemales = ["I_F1-1_F", "I_F16-2_F", "I_F9-1_F", "I_G13-2_F", "I_GX_F", "I_S4-19_M", "I_LA_F", "I_LY-1_F", "I_S3_F" ,"I_S2_F","I_GY-2_F"],\
	species ="collared")


OC =  SBEpop(
	popcode = "OC",
	vcf_snps =  "/home/ludovic/nobackup_ludo/SBE/source_files/gt_OC.rm.vcf.gz",
	vcf_allsites = "/home/ludovic/nobackup_ludo/SBE/source_files/OC.gatk.allsites.rm.vcf.gz",
	vcf_snps18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/gt_OC_forfst.rm.vcf.bgz",
	males9 = ["OC_1_M", "OC_2_M", "OC_3_M", "OC_4_M", "OC_5_M", "OC_HB1_M", "OC_HB2_M", "OC_HB3_M", "OC_HB5_M" ],\
	females9 = ["OC_HB6_F", "OC_HB7_F", "OC_HB8_F", "OC_HB9_F", "OC_HB10_F", "OC_10_F", "OC_6_F", "OC_7_F", "OC_8_F"],
	allmales = ["OC_1_M", "OC_2_M", "OC_3_M", "OC_4_M", "OC_5_M", "OC_HB1_M", "OC_HB2_M", "OC_HB3_M", "OC_HB5_M" ,"OC_HB4_M"],\
	allfemales = ["OC_HB6_F", "OC_HB7_F", "OC_HB8_F", "OC_HB9_F", "OC_HB10_F", "OC_10_F", "OC_6_F", "OC_7_F", "OC_8_F"],	
	species ="collared")


H =  SBEpop(
	popcode = "H",
	vcf_snps =  "/home/ludovic/nobackup_ludo/SBE/source_files/gt_H.rm.vcf.gz",
	vcf_allsites = "/home/ludovic/nobackup_ludo/SBE/source_files/H.gatk.allsites.rm.vcf.gz",
	vcf_snps18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/gt_H_forfst.rm.vcf.bgz",
	males9 = ["H_117_M", "H_12_M", "H_367_M", "H_379_M", "H_452_M", "H_454_M", "H_94_M", "H_76_M", "H_86_M" ],\
	females9 = ["H_118_F", "H_126_F", "H_31_F", "H_354_F", "H_377_F", "H_392_F", "H_43_F", "H_78_F", "H_52_F" ],\
	allmales = ["H_117_M", "H_12_M", "H_367_M", "H_379_M", "H_452_M", "H_454_M", "H_94_M", "H_76_M", "H_86_M" ,"H_53_M"],\
	allfemales = ["H_118_F", "H_126_F", "H_31_F", "H_354_F", "H_377_F", "H_392_F", "H_43_F", "H_78_F", "H_52_F","H_473_F" ],\
	species ="collared")
