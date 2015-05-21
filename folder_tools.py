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
		If exceptions==True : an exception is raised if any  job has not properly been terminated

#TESTS
>>> check_logs("test_files",exceptions = False)
Found 1 file(s) and  1 file(s) that are unfinished or with error
error_files:
test_files/test.log
>>> check_logs("test_files",".log2" ,exceptions = False)
Found 1 file(s) and  0 file(s) that are unfinished or with error
		"""
	logfiles = [folder+"/"+filename for filename in os.listdir(folder) if filename.endswith(pattern) ]
	n_logs_files = 0
	files_fail =[] 
	for filename in logfiles:
		n_logs_files +=1
		lastline = sh.tail("-n",1,filename)
		if not "SLURM: end" in lastline: 
			files_fail.append(filename)
			#print "ERROR: "+filename+ " \t endswith : " + str(lastline)
	print "Found "+ str(n_logs_files)+ " file(s) and  "+ str(len(files_fail)) +" file(s) that are unfinished or with error"
	if  len(files_fail) != 0:
		print "error_files:\n"+"\n".join(files_fail)

