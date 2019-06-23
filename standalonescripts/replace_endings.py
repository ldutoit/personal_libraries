#!/usr/bin/env python
##EXAMPLE USAGE: python replace_endings.py MYFILE
import sys
output=open(sys.argv[1]+"_withendings.fasta","w")
with open(sys.argv[1]) as f:
	for line in f:
		if line.startswith(">"):
			output.write(line)
		else:
			output.write(line.strip()+"*\n")

output.close()

