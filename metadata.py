#!/usr/bin/env python2
# Filename:metadata_SBE.py


##This library use different metadata I made used of throughout my PhD working with various flycatcher datasets

class SBEpop(object):
	''' a class that has file location and sex information for the project about Sex biased gene expression. 
	males9 and females 9 are the 9 individuals with the highest coverage, they are used for FST calculations while all
	individuals are used for pi
	self.ind18 is the sum of males9 and feamles 9.  '''
	def __init__(self,popcode,vcf_snps,vcf_allsites,vcf_snps18,vcf_allsites18,males9,females9,allmales,allfemales,species):
		self.popcode = popcode
		self.vcf_snps = vcf_snps
		self.vcf_allsites = vcf_allsites
		self.vcf_snps18 = vcf_snps18
		self.vcf_allsites18 = vcf_allsites18
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
	vcf_allsites18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/CZC.gatk.allsites_forfst.rm.vcf.gz",
	males9 = ["CZC_304236_M", "CZC_382761_M", "CZC_419135_M", "CZC_83359_M", "CZC_48837_M", "CZC_82259_M", "CZC_82355_M", "CZC_82356_M", "CZC_83353_M" ],\
	females9 = ["CZC_119601_F", "CZC_382648_F", "CZC_419251_F", "CZC_91360_F", "CZC_73077_F", "CZC_82317_F", "CZC_82348_F", "CZC_82352_F", "CZC_91319_F" ],\
	allmales = ["CZC_304236_M", "CZC_382761_M", "CZC_419135_M", "CZC_83359_M", "CZC_48837_M", "CZC_82259_M", "CZC_82355_M", "CZC_82356_M", "CZC_83353_M" ,"CZC_461398_M"],\
	allfemales =["CZC_119601_F", "CZC_382648_F", "CZC_419251_F", "CZC_91360_F", "CZC_73077_F", "CZC_82317_F", "CZC_82348_F", "CZC_82352_F", "CZC_91319_F" , "CZC_419493_F"],\
	species ="collared")


I =  SBEpop(
	popcode = "I",
	vcf_snps =  "/proj/b2010010/repos/variation/200genomes/vcf/population_genotypes/R1I1_unfiltered/gt_I.vcf.gz",
	vcf_allsites = "/proj/b2010010/repos/variation/200genomes/vcf/allsites/I.gatk.allsites.rm.vcf.gz",
	vcf_snps18 = "",
	vcf_allsites18 = "",
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
	vcf_allsites18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/OC.gatk.allsites_forfst.rm.vcf.gz",
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
	vcf_allsites18 = "/home/ludovic/nobackup_ludo/SBE/vcfs/H.gatk.allsites_forfst.rm.vcf.gz",
	males9 = ["H_117_M", "H_12_M", "H_367_M", "H_379_M", "H_452_M", "H_454_M", "H_94_M", "H_76_M", "H_86_M" ],\
	females9 = ["H_118_F", "H_126_F", "H_31_F", "H_354_F", "H_377_F", "H_392_F", "H_43_F", "H_78_F", "H_52_F" ],\
	allmales = ["H_117_M", "H_12_M", "H_367_M", "H_379_M", "H_452_M", "H_454_M", "H_94_M", "H_76_M", "H_86_M" ,"H_53_M"],\
	allfemales = ["H_118_F", "H_126_F", "H_31_F", "H_354_F", "H_377_F", "H_392_F", "H_43_F", "H_78_F", "H_52_F","H_473_F" ],\
	species ="collared")




oland96 = SBEpop(
	popcode = "OC96",
	vcf_snps =  "source_files/gt_Ficedula_200plus90.vcf.gz",
	vcf_allsites = "vcfs/oland96_allsites.gz",
	vcf_snps18 = "NA",
	vcf_allsites18 = "NA",
	males9 = "NA",\
	females9 = "NA",\
	allmales =  ["B475143","B476196","u001", "u002", "u003", "u004", "u005", "u007", "u010", "u011", "u012", "u013", "u014", "u016", "u017", "u018", "u019", "u020", "u021", "u022", "u023", "u024", "u025", "u026", "u027", "u028", "u029", "u030", "u031", "u032", "u033", "u034", "u036", "u037", "u038", "u040", "u041", "u042", "u043", "u101", "u102", "u103", "u104", "u105", "u106", "u107", "u108", "u109", "u110", "u111", "u112", "u113", "u114", "u115", "u116", "u117", "u118", "u119", "u120", "u121", "u122", "u123", "u125", "u126", "u128", "u129", "u130", "u131", "u132", "u133", "u134", "u135", "u136", "u137", "u138", "u139", "u140", "u141", "u142", "u143", "u144", "u145", "u146", "u147", "u148", "u149"] +  ["OC_1_M", "OC_2_M", "OC_3_M", "OC_4_M", "OC_5_M", "OC_HB1_M", "OC_HB2_M", "OC_HB3_M", "OC_HB5_M", "OC_HB6_F"] ,\
	allfemales = [], \
	species ="collared") 
gotland96 =  SBEpop(
	popcode = "GO96",
	vcf_snps =  "vcfs/gotland96.SNP.vcf.gz",
	vcf_allsites = "vcfs/gotland97_allsites.gz",
	vcf_snps18 = "",
	vcf_allsites18 = "",
	males9 = [],\
	females9 = [],\
	allmales =  [ "15M153", "15M155", "15M158", "15M160", "15M161", "15M162", "15M163", "15M201", "15M202", "15M203", "15M204", "15M207", "15M468", "15M469", "15M475", "15M477", "15M49", "15M537", "15M568", "15M571", "15M573", "15M589", "15M684", "15M724", "93M25", "93M27", "93M28", "93M29", "93M36", "93M38", "93M39", "93M40", "93M41", "93M46", "93M53", "93M55", "93M58", "93M71", "93M72", "93M73", "93M78", "93M79", "93M80", "93M81", "93M83", "93M84", "93M86", "93M91"],\
	allfemales = [ "15F129", "15F130", "15F131", "15F135", "15F142", "15F143", "15F145", "15F149", "15F151", "15F17", "15F18", "15F21", "15F22", "15F23", "15F24", "15F25", "15F29", "15F447", "15F448", "15F450", "15F453", "15F457", "15F459", "15F460", "93F24", "93F26", "93F30", "93F31", "93F32", "93F34", "93F35", "93F42", "93F44", "93F45", "93F47", "93F54", "93F56", "93F59", "93F74", "93F75", "93F77", "93F82", "93F88", "93F89", "93F90", "93F92", "93F93", "93F94"],\
	species ="collared") 
