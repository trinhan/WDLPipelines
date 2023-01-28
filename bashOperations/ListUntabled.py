'''
Script to rename name paths in a Terra tsv file to a new directory
'''
import os
from csv import reader
import argparse
import re
from re import match
import pandas as pd
import numpy as np
import subprocess

def return_dict(fields, values):
	''' Function takes a list of fields of interest e.g. clinRep and matches to the appropriate
		header in the tsv file. Returns a dictionary of string: matched headers
	'''
	test=[list(filter(lambda v: match(x, v), values)) for x in fields]
	d=dict(zip(fields, test))
	return d



def process_fle(inTsv,mapper, dest, outTsv):
	''' Function an input tsv and reads each line. The first line is used to create a mapper
		dictionary of which fields to change. 
		Each line is read in. If the column name matches a string in the dict, change the file path to
		gs://dest/mapperstring/file. This is then written to outTsv file
	'''
	orig_list=[]
	new_list=[]
	with open (outTsv, 'w') as writer:
		with open(inTsv) as f:
			for i, data in enumerate(reader(f, delimiter="\t")):
				if i==0:
					''' header read in here and dictionary created
					'''
					g= open(mapper, "r")
					fields =  g.read().splitlines()
					d= return_dict(fields, data)
					headerlist=data
					writer.write('\t'.join(data)+'\n')
				if i>0:
					''' Each line read in individually, change the file paths
					'''
					for x in d:
						for y in d[x]:
							nx=headerlist.index(y)
							var=data[nx].split("/")
							keep=var[len(var)-1]
							''' Only rename the path if the cell is not empty
							'''
							if keep!="" :
								newVar="gs://"+dest+"/"+x+"/"+keep
								orig_list.append(data[nx])
								data[nx]=newVar
								new_list.append(newVar)
					writer.write('\t'.join(data)+'\n')
	writer.close()
	f.close()
	ChangeList=zip(orig_list, new_list)
	with open('oldList.txt', 'w') as outList:
		outList.write('\n'.join(orig_list))
	outList.close()
	with open('newList.txt', 'w') as outList:
		outList.write('\n'.join(new_list))
	outList.close()

def greppattern(ar,pattern):
    indices = [i for i,x in enumerate(ar) if re.match(pattern,x)]
    print(indices)
    return ar.take(indices)


def collapseList(inTsv,searchTerm):
	# initializing datelist
	df=pd.read_csv(inTsv, sep="\t")
	df2=df.values.flatten()
	df3=pd.Series(df2)
	mask = df3.str.contains(searchTerm)
	outlist=df3[mask==True]
	return (outlist)


def searchDirectory(bucket,searchTerm, SafeList, outputfile):
	# issue a gsutil command to search
	command = (f"gsutil ls gs://{bucket}/**/*.{searchTerm}")
	delList=os.popen(command).read().strip()
	delList = delList.split("\n")
	# compare to the main list:
	main_list=list(set(delList) - set(SafeList))
	print('The following samples not on safe list:')
	print(main_list)
	file = open(outputfile,'w')
	for item in main_list:
		file.write(item+"\n")
	file.close()
	print('Items are written to file! Finished!')



def parseOptions():
# The main argument parser
	parser = argparse.ArgumentParser(description="Script for changing file locations from a Terra workspace")
	parser.add_argument('-i', '--inputtsv', help='input tsv file (download samples table in Terra)')
	parser.add_argument('-d', '--bucket', help="directory to list files in")
	parser.add_argument('-s', '--searchTerm', help="term to search for e.g. 'bam'") 
	parser.add_argument('-o', '--outputtsv', help="Name of output file") 
	args = parser.parse_args()
	return args

def main():
	args = parseOptions()
	if(args):
		SafeList=collapseList(args.inputtsv,args.searchTerm)
		PrintList=searchDirectory(args.bucket,args.searchTerm,SafeList, args.outputtsv)
	else:
		parser.print_help()

if __name__ == "__main__":
	main()
