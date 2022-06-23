'''
Script to rename name paths in a Terra tsv file to a new directory
'''
import os
from csv import reader
import argparse
import re
from re import match

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
							orig_list.append(data[nx])
							var=data[nx].split("/")
							keep=var[len(var)-1]
							''' Only rename the path if the cell is not empty
							'''
							if keep!="" :
								newVar="gs://"+dest+"/"+x+"/"+keep
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



def parseOptions():
# The main argument parser
	parser = argparse.ArgumentParser(description="Script for changing file locations from a Terra workspace")
	parser.add_argument('-i', '--inputtsv', help='input tsv file (download samples table in Terra)')
	parser.add_argument('-m', '--mapper', help="mapper: list of strings to search for. Matching samples will be moved to dest/mapper[x]")
	parser.add_argument('-d', '--destination', help="destination google bucket WITHOUT gs:// ")
	parser.add_argument('-o', '--outputtsv', help="output tsv name (for reupload to Terra workspace") 
	args = parser.parse_args()
	return args

def main():
	args = parseOptions()
	if(args):
		process_fle(args.inputtsv, args.mapper, args.destination, args.outputtsv)
	else:
		parser.print_help()

if __name__ == "__main__":
	main()
