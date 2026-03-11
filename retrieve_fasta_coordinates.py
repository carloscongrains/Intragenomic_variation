#!/usr/bin/env python

"""
Retrieve fasta sequences from a list of ids and coordinates and a multifasta file (database).
If first coordinate greater than second coordinate in the input_list file, the program will generate a reverse complement sequence.
"""

import argparse,os,sys
from Bio import SeqIO
from Bio.Seq import Seq
import copy

class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_list', help='Input ids and coordinates in a tab separated file.')
parser.add_argument('--input_fasta', help='Path of the fasta file. One ID (first word of the fasta header) per line.')
parser.add_argument('--output_fasta', help='Path of the output file')
parser.add_argument('--mode', help='Choose between "fast" or "slow". If the input fasta file is huge use slow to prevent run out of RAM.')

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()

if args.input_list:
	input_list = args.input_list

if args.input_fasta:
	input_fasta = args.input_fasta

if args.output_fasta:
	output_fasta = args.output_fasta
	
if args.mode:
	mode = args.mode

assert mode == "fast" or mode == "slow", "Mode must be: fast or slow"
assert os.stat(input_list).st_size >= 0.0, "Your input list is empty"

"""
Main
"""

#Preparing the list containing the ids.

id_list_coordnates = []
with open(input_list) as input_list_handle:
	lines = input_list_handle.readlines()

for info in lines:
	if not info.isspace():
		info = info.rstrip()
		info_list = info.split("\t")
		coord_1 = int(info_list[1])
		coord_2 = int(info_list[2])
		id_list_coordnates.append([info_list[0],coord_1,coord_2])

#Processing and saving the results
if mode == "fast":
	record_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
	with open(output_fasta, "w") as output_handle:
		for ID_cord in id_list_coordnates:
			if ID_cord[0] in record_dict:
				record_info = record_dict[ID_cord[0]]
				record = copy.deepcopy(record_info)
				seq = str(record.seq)
				if ID_cord[2]>ID_cord[1]:
					if ID_cord[1] > 0 and ID_cord[2] <= len(seq):
						new_seq = seq[ID_cord[1]-1:ID_cord[2]]
						record.seq = Seq(new_seq)
						new_id = record.id + "_" + str(ID_cord[1]) + "-" + str(ID_cord[2])
						record.id = new_id
						record.description = ""
						record.name = ""
						output_handle.write(record.format("fasta"))
					else:
						print(ID_cord[0] + " slice is not in the record (length = " + str(len(seq)) + "). Check the coordinates " + str(ID_cord[1]) + " and " + str(ID_cord[2]) + ". " + input_fasta)
				elif ID_cord[2]<ID_cord[1]:
					if ID_cord[2] > 0 and ID_cord[1] <= len(seq):
						new_seq = seq[ID_cord[2]-1:ID_cord[1]]
						record.seq = Seq(new_seq).reverse_complement()
						new_id = record.id + "_revcomp_" + str(ID_cord[1]) + "-" + str(ID_cord[2])
						record.id = new_id
						record.description = ""
						record.name = ""
						output_handle.write(record.format("fasta"))
					else:
						print(ID_cord[0] + " slice is not in the record (length = " + str(len(seq)) + "). Check the coordinates " + str(ID_cord[1]) + " and " + str(ID_cord[2]) + ". " + input_fasta)
			else:
				print(ID_cord[0] + " is not in " + input_fasta)

elif mode == "slow":
	record_dict = SeqIO.index(input_fasta, "fasta")
	with open(output_fasta, "w") as output_handle:
		for ID_cord in id_list_coordnates:
			if ID_cord[0] in record_dict:
				record = record_dict[ID_cord[0]]
				seq = str(record.seq)
				if ID_cord[2]>ID_cord[1]:
					if ID_cord[1] > 0 and ID_cord[2] <= len(seq):
						new_seq = seq[ID_cord[1]-1:ID_cord[2]]
						record.seq = Seq(new_seq)
						new_id = record.id + "_" + str(ID_cord[1]) + "-" + str(ID_cord[2])
						record.id = new_id
						record.description = ""
						record.name = ""
						output_handle.write(record.format("fasta"))
					else:
						print(ID_cord[0] + " slice is not in the record (length = " + str(len(seq)) + "). Check the coordinates " + str(ID_cord[1]) + " and " + str(ID_cord[2]) + ". " + input_fasta)
				elif ID_cord[2]<ID_cord[1]:
					if ID_cord[2] > 0 and ID_cord[1] <= len(seq):
						new_seq = seq[ID_cord[2]-1:ID_cord[1]]
						record.seq = Seq(new_seq).reverse_complement()
						new_id = record.id + "_revcomp_" + str(ID_cord[1]) + "-" + str(ID_cord[2])
						record.id = new_id
						record.description = ""
						record.name = ""
						output_handle.write(record.format("fasta"))
					else:
						print(ID_cord[0] + " slice is not in the record (length = " + str(len(seq)) + "). Check the coordinates " + str(ID_cord[1]) + " and " + str(ID_cord[2]) + ". " + input_fasta)
			else:
				print(ID_cord[0] + " is not in " + input_fasta)
