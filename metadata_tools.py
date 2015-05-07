#!/usr/bin/env python2
# Filename: metadata_tools.py
#small class and objects use to manipulate rapidly inds and populations
import pandas as pd
import matplotlib.pyplot as plt
import vcf

path_to_metadata="/home/ludovic/private/ind_metadata.txt"
inds_df=pd.read_csv(path_to_metadata,sep="\t",na_values="NA")# data frame of the metadata file

class Ind(object):
	""" a class to define ind, contain information in the metadata file
	""" 
	def __init__ (self,name):
		self.NAME=str(inds_df["NAME"][inds_df["NAME"]== name]).split()[1]
		self.SPECIES=str(inds_df["SPECIES"][inds_df["NAME"]== name]).split()[1]
		self.POP=str(inds_df["POP"][inds_df["NAME"]== name]).split()[1]
		self.SEX=str(inds_df["SEX"][inds_df["NAME"]== name]).split()[1]
		self.POP_SPE=str(inds_df["POP_SPE"][inds_df["NAME"]== name]).split()[1]
		self.LATIN_SPE=str(inds_df["LATIN_SPE"][inds_df["NAME"]== name]).split()[1]
		###cov to deal with NA
		MEDIAN=str(inds_df["COV_MEDIAN"][inds_df["NAME"]== name]).split()[1]
		MEAN=str(inds_df["COV_MEAN"][inds_df["NAME"]== name]).split()[1]
		SD=str(inds_df["COV_SD"][inds_df["NAME"]== name]).split()[1]
		###vcf file
		self.VCF_FOLDER=str(inds_df["VCF_FOLDER"][inds_df["NAME"]== name]).split()[1]
		self.VCF=str(inds_df["VCF"][inds_df["NAME"]== name]).split()[1]
		if MEDIAN != "NaN": 
			self.COV_MEDIAN=int(MEDIAN)
			self.COV_MEAN=int(MEAN)
			self.COV_SD=int(SD)
		else:
			self.COV_MEDIAN=self.COV_MEAN=self.COV_SD=None
	def __repr__(self):
		toreturn="    ".join([str(self.NAME),str(self.SPECIES),str(self.POP),str(self.SEX)])
		return repr(toreturn)



#dict of ind where keys are names and values are Ind object
dict_inds={}
for item in inds_df["NAME"]:
	dict_inds[item]=Ind(item)


##colors dictionnary for plots
pop_colors={
'italy':'red',
'czech':'green',
'indonesia': 'gold',
'serbia': 'blue',
'morocco': 'orange',
'sweden': 'gold',
'oland': 'black',
'hungary': 'yellow',
'spain': 'purple',
'bulgaria': 'aqua'}





pop_spe_colors={
 'par_swe':"indigo",
 'hyp_ind' :"purple",
 'hyb_ola':"pink",
 'sem_bul':"brown",
 'atl_mor':"grey",

 'col_hun':"darkgreen",
 'col_ser':"lightblue",
 'col_ola':"cyan",
 'col_ita':"green",
 'col_cze':"blue",

 'pie_cze':"orange",
 'pie_ola':"gold",
 'pie_spa':"red",
 'pie_swe':"darkred",
 }

# s

