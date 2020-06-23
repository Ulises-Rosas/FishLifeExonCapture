#!/usr/bin/env python
import io
import re
import argparse

from Bio.SeqIO import to_dict, parse, write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def getOpts(): 
	parser = argparse.ArgumentParser(description="""
										Filters the exonerate output for sequences on the
										correct strand, and renames the sequence simply 
										as >Taxon""")
	parser.add_argument('-f', '--fasta' ,
						dest = 'fasta' , 
						type = str ,
						default= None ,
						required= True, 
						help = 'Fasta alignment to prune')
	parser.add_argument('-o', '--output', 
						dest = 'output',
						type = str, 
						default = None, 
						required = True,
						help = 'Name of output file')
	parser.add_argument('-t', '--taxon',
						dest = 'taxon',
						type = str, 
						default = None, 
						required = True,
						help = 'Name of taxon')
	parser.add_argument('-l', '--flanks',
						dest = 'flanks',
						type = str,
						default = None,
						required = False,
						help = 'Complete assembled sequences')
	args, unknown = parser.parse_known_args()
	return args

def matchExonerateCdhit(fasta, flanks, taxon, output):
	""" 
	my %opts = { fasta  => "*cdhit.exonerate",
	    		 flanks => "*.cdhit",
				 taxon  => $core }
	"""

	# read in sequences
	input = open(fasta, 'r').readlines()
	input = [i for i in input if not re.findall('(Command|Hostname|exonerate)', i)]
	input = io.StringIO( "".join(input) )
	records = list(parse(input, "fasta"))

	# filter just the sequences in the correct orientation 
	# (known from the reference sequence used with exonerate)
	oriented = []
		
	for record in records:
		coordinates = record.description.split("\t")[1]
		start = int(coordinates.split("-")[0])
		end = int(coordinates.split("-")[1])
		if end > start:
			oriented.append(record)	

	WSdict  = to_dict(parse(flanks,"fasta"))
	genomic = []
		
	for sequence in oriented:
		genomicSeq = WSdict[sequence.id]
		if sequence.seq in genomicSeq.seq:
			newRecord = SeqRecord(genomicSeq.seq, id=taxon, description='')
			genomic.append(newRecord)
		else:
			genomicSeqRev = genomicSeq.seq.reverse_complement()	
			newRecord = SeqRecord(genomicSeqRev, id=taxon, description='')
			genomic.append(newRecord)
				
	write(genomic, output, "fasta")

def main():

	args = getOpts()

	filterFlanks(fasta  = args.fasta,
				 output = args.output,
				 taxon  = args.taxon,
				 flanks = args.flanks)

if __name__ == '__main__':
	main()
