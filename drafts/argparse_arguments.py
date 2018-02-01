#!/usr/bin/env python
#import modules
import argparse

#parser
parser = argparse.ArgumentParser() # add the parse
parser.add_argument("square", help="display a square of a given number", type=int) # example of a non string positional argument, it is good to specify type
parser.add_argument("--verbose", help="increase output verbosity", action="store_true") #An optional argument (or option) is (by default) given None as a value when its not being used. The keyword "action" is being given the value "store_true" which means that if the option is specifed, then assign the value "True" to args.verbose
args = parser.parse_args()	# returns data from the options specified (infile)




## program
print(args.square**2)	

if args.verbose:
    print("verbosity turned on")