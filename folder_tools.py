#!/usr/bin/python
# Filename: folder_tools.py
import os,sh

# class File(object):
# 	def __init__ (self,path):
# 		self.name=os.path.basename(path)
# 		self.abs=os.path.abspath(name)
# 		self.ext=os.path.splitext(name)[1]
# 		self.root_abs=os.path.splitext(self.abs)[0]
# 		self.root_rel=os.path.splitext(name)[0]
# 		self.folder=
# 		self.exist=
# 		self.size=


def check_logs(folder,pattern = ".log",exceptions=True):
	"""Check all the logfiles in a folder designated by a pattern for the presence of 'SLURM: end'.
		If exceptions==True : an exception is raised if any  job has not properly been terminated"""
	logfiles = [folder+"/"+filename for filename in os.listdir(folder) if pattern in filename ]
	n_logs_files = 0
	n_logs_fail = 0
	files_fail =[] 
	for filename in logfiles:
		n_logs_files +=1
		lastline = sh.tail("-n",1,filename)
		if not "SLURM: end" in lastline: 
			n_logs_fail +=1
			files_fail.append(filename)
			print "ERROR: "+filename+ " \n"+ str(lastline) +"\n"
	print "checked "+ str(n_logs_files)+ " "+ str(n_logs_fail) +" unfinished/error jobs"
	if exceptions and n_logs_fail != 0:
		print "error_files:\n"+"\n".join(files_fail)
		raise Exception

