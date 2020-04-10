#!/usr/bin/env python

import argparse
from Bio import SeqIO

def getOpt():
	parser = argparse.ArgumentParser(description="Get the longest")
	parser.add_argument('-f', '--fasta' , dest = 'fasta' , type = str , default= None , required= True, help = 'Fasta alignment to prune')
	parser.add_argument('-o', '--output', dest = 'output', type = str, default = None, required = True, help = 'Name of output file')
	args, unknown = parser.parse_known_args()
	return args

def writeLong(fasta, output):
	records = list(SeqIO.parse(fasta, "fasta"))
	records.sort(key=lambda r: -len(r))
	SeqIO.write(records[0], output, "fasta")

def main():
	args = getOpt()
	writeLong(args.fasta, args.output)

if __name__ == '__main__':
	main()
