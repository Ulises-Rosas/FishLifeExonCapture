#!/usr/bin/env python

import io
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def getOpts():
	parser = argparse.ArgumentParser(description="""

		Filters the exonerate output for sequences on the correct strand, and renames the sequence simply as >Taxon

		""")
	parser.add_argument('-f', '--fasta',
						dest = 'fasta' ,
						type = str ,
						default= None ,
						required= True,
						help = 'Fasta alignment to prune')
	parser.add_argument('-o','--output',
						dest = 'output',
					    type = str, 
					    default = None,
					    required = True,
					    help = 'Name of output file')
	parser.add_argument('-t',
					    '--taxon',
					     dest = 'taxon',
					     type = str, 
					     default = None,
					     required = True, 
					     help = 'Name of taxon')

	args, unknown = parser.parse_known_args()
	return args
	

def getoriented(fasta, taxon, output):

	# read in sequences
	input = open(fasta, 'r').readlines()
	input = [i for i in input if not re.findall('(Command|Hostname|exonerate)', i)]
	input = io.StringIO( "".join(input) )
	records = list(SeqIO.parse(input, "fasta"))

	# filter just the sequences in the correct orientation (known from the reference sequence used with exonerate)
	oriented = []

	for record in records:
		coordinates = record.description.split("\t")[1]
		start = int(coordinates.split("-")[0])
		end = int(coordinates.split("-")[1])
		if end > start:
			oriented.append(record)	

	filtered = []

	for record in oriented:
		filteredSeq = SeqRecord(record.seq, id=taxon, description='')
		filtered.append(filteredSeq)

	SeqIO.write(filtered, output, "fasta")

def main():

	args = getOpts()
	getoriented(args.fasta, args.taxon, args.output)

if __name__ == '__main__':
	main()
