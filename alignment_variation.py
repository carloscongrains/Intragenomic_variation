#!/usr/bin/env python

"""
This script calculates: 
1. number of sequences
2. number of sequences without ambiguities
3. number of sequences with ambiguities
4. alignment length	
5. Number of variable sites including nucleotides and indels
6. Number of variable sites including only nucleotides
7. Number of parsimony informative sites	
8. Number of variable sites including only indels
9. Mean of pairwise distance including gaps (gapweight=0.5)	calculated using distmat tool (https://www.bioinformatics.nl/cgi-bin/emboss/distmat)
10. Maximum pairwise distance including gaps (gapweight=0.5)	calculated using distmat tool (https://www.bioinformatics.nl/cgi-bin/emboss/distmat)
"""

import argparse,os,sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from collections import Counter
import subprocess
import numpy
import csv


class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_fasta', help='Path of the fasta dir. This directory must contain subdiretories of molecular markers. Each subdirectory must have fasta files with names: samplesID.fas.')
#parser.add_argument('--output', help='A tab delimited file to store the results. It contains the following fields: markerID, sampleID, alignment length, number of matches, number of mismatches, number of single nucleotide mismatches, number of gap mismatches.')


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()

if args.input_fasta:
	input_fasta = args.input_fasta

#if args.output:
#	output = args.output

	
"""
FUNCTIONS
"""

#Convert relative to full path
def convert_to_full_path(path):
    #Get working directory
    absolute_path = os.path.abspath(path)
    return absolute_path

#Alignment function if needed
def alignment(fasta_path,output_path):

	output_dir = os.path.dirname(output_path)
	infile = fasta_path
	outfile = output_path
	with open(outfile,"wb") as out_handle:
		cmd = ["mafft", "--localpair","--maxiterate","1000","--thread","1","--quiet","--preservecase", infile]
		p = subprocess.Popen(cmd, stdout=out_handle,stderr=subprocess.PIPE, cwd=output_dir)
		p.wait()
		(stdout, stderr) = p.communicate()
	#	except subprocess.calledProcessError as err:
	#		log.info("alignment failed - {} - {}".format(err.stderr, fasta_path))
	#		return
	return output_path

def trimming(input_align):
	align = AlignIO.read(input_align, "fasta")
	aln_len = align.get_alignment_length()

	seq_slice = []
	for pos in range(aln_len):
		subseq = align[:, pos]
		if "-" in subseq or "N" in subseq:
			continue
		else:
			seq_slice.append(pos)
			break

	#Check the alignemnt backwards
	for pos in range(aln_len-1,-1,-1):
		subseq = align[:, pos]
		if "-" in subseq or "N" in subseq:
			continue
		else:
			pos_end = pos + 1
			seq_slice.append(pos_end)
			break

	trimming_aln = align[:,seq_slice[0]:seq_slice[1]]
	number_seq = len(align)
	return trimming_aln,number_seq

#Belongs to wraper_cleaning_alingment function
#Get the positions with gap-only or N-only
def uniformative_sites(input_Bioalign):
	#input_Bioalign = AlignIO.read(input_align, "fasta")
	aln_len = input_Bioalign.get_alignment_length()
	position_to_remove = []
	seq_slice = []
	for pos in range(aln_len):
		subseq = input_Bioalign[:, pos]
		subseq = subseq.upper()
		subseq_list = list(subseq)
		subseq_list_uniq = list(set(subseq_list))

		if ["-"] == subseq_list_uniq or ["N"] == subseq_list_uniq:
			position_to_remove.append(pos)

	return position_to_remove

#Belongs to wraper_cleaning_alingment function
#Remove sites with gap-only or N-only
def edit_sequences(input_Bioalign,position_to_remove):
	position_to_remove.sort(reverse=True)
	list_sequence = []
	for record in input_Bioalign:
		sequence = str(record.seq).upper()
		sequence_list = list(sequence)
		for position in position_to_remove:
			del(sequence_list[position])
		sequence = "".join(sequence_list)
		list_sequence.append([record.id,sequence])
	length_aln = len(list_sequence[0][1])
	return len(position_to_remove),length_aln, list_sequence

#Belongs to wraper_cleaning_alingment function
#Save the results
def save_output(list_sequence,output):
	with open(output,"w") as out_handle:
		for record in list_sequence:
			out_handle.write(">" + record[0] + "\n" + record[1] + "\n")
	return

#Function to remove sites with gap-only or N-only
def wraper_cleaning_alingment(input_Bioalign,output):
	position_to_remove = uniformative_sites(input_Bioalign)
	len_list_sequence,length_aln, list_sequence = edit_sequences(input_Bioalign,position_to_remove)
	save_output(list_sequence,output)
	return length_aln

#Fasta alignment to Bioseq alignment
def parse_fasta(input_align):
	align = AlignIO.read(input_align, "fasta")
	number_seq = len(align)
	return align,number_seq

#Remove sequences IUPAC code for 3 or more possible nucleotides 
def check_amb(alignment):
	#"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","T","C","G"]
	amb_iupac = ["B","D","H","V","N","M","R","W","S","Y","K"]
	number_seq = len(alignment)
	
	#Create new alignment
	#new_align = MultipleSeqAlignment([])
	multiple_alignment_list = []
	number_seq_with_ambiguities = 0
	control = 0
	for idx in range(number_seq):
		record = alignment[idx]
		for amb_base in amb_iupac:
			if amb_base in str(record.seq):
				number_seq_with_ambiguities = number_seq_with_ambiguities + 1
				control = 1
				break
		if control == 0:
			record_aln = SeqRecord(record.seq, id=record.id, description="")
			multiple_alignment_list.append(record_aln)
		control = 0

	number_seq_without_ambiguities = len(multiple_alignment_list)
	new_align = MultipleSeqAlignment(multiple_alignment_list)
	
	return new_align,number_seq_without_ambiguities,number_seq_with_ambiguities

#Save alignment into a file
def aln_to_fasta_file(align,output_file):
	with open(output_file, 'w') as output_handle:
		for record in align:
			output_handle.write(">" + record.id + "\n" + str(record.seq) + "\n")
	return

#Analyze the alignment to get number of variable
def seq_comparison(alignment):

	length_aln = alignment.get_alignment_length()
	
	counter_variable = 0
	counter_variable_nuc = 0
	counter_variable_nuc_parsimony_inf = 0
	counter_variable_indel = 0

	for pos in range(length_aln):
		seq_pos = alignment[:, pos]
		seq_pos_upper = seq_pos.upper()
		seq_pos_freq = Counter(seq_pos_upper)
		if len(seq_pos_freq.keys()) == 1:
			continue
		elif len(seq_pos_freq.keys()) == 2:
			counter_variable = counter_variable + 1
			if "-" in seq_pos_freq.keys():
				counter_variable_indel = counter_variable_indel + 1
			else:
				counter_variable_nuc = counter_variable_nuc + 1
				if parsimony_informative_site(seq_pos_freq) == "yes":
					counter_variable_nuc_parsimony_inf = counter_variable_nuc_parsimony_inf + 1

		elif len(seq_pos_freq.keys()) >= 3:
			counter_variable = counter_variable + 1
			counter_variable_nuc = counter_variable_nuc + 1
			if parsimony_informative_site(seq_pos_freq) == "yes":
				counter_variable_nuc_parsimony_inf = counter_variable_nuc_parsimony_inf + 1
			if "-" in seq_pos_freq.keys():
				counter_variable_indel = counter_variable_indel + 1

	return counter_variable, counter_variable_nuc, counter_variable_nuc_parsimony_inf, counter_variable_indel

#Get parsimony-informative sites
def parsimony_informative_site(seq_pos_freq):
	counter = 0
	for base in seq_pos_freq.keys():
		if base in ["A","T","C","G"]:
			if seq_pos_freq[base] >= 2:
				counter = counter + 1
				if counter == 2:
					return "yes"
	return "no"


#Main function to get the consensus

def distmat(path_alignment_fasta,input_DIR):
	"""Convert the distmat result to a distance matrix in phylip format
	"""

	working_dir = input_DIR

	#Create the path to trimal dna dmat
	path_alignment_fasta_dmat = os.path.splitext(path_alignment_fasta)[0] + ".dist"

	#Run distmat
	cmd = ["distmat","-nucmethod","0","-sequence",path_alignment_fasta,"-outfile",path_alignment_fasta_dmat,"-gapweight","0.5"] 
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, cwd= working_dir)
	out = p.communicate()

	return path_alignment_fasta_dmat


def distmat_to_np_array(path_alignment_fasta_dmat):
	"""Parse distmat output to a numpy array
	"""

	#Open the dist mat file
	with open(path_alignment_fasta_dmat,"r") as dist_handle:
		lines_dist_matrix = dist_handle.readlines()

	#Create the list to put the results
	distance_matrix_final=[]
	list_ids = []

	#Populate the list of distance matrix and list of ids
	for index_dist in range(8,len(lines_dist_matrix)):
		list_dist = lines_dist_matrix[index_dist].split("\t")[1:-2]
		list_ids.append(lines_dist_matrix[index_dist].split("\t")[-1].strip().split(" ")[0])
		for index_pair_dist in range(len(list_dist)):
			if list_dist[index_pair_dist] == '':
				list_dist[index_pair_dist] = float(0.00)
			else:
				list_dist[index_pair_dist] = float(list_dist[index_pair_dist])
		distance_matrix_final.append(list_dist)

	#Convert into a numpy array
	distance_matrix_array = numpy.array(distance_matrix_final)
	num_col = numpy.size(distance_matrix_array,1)
	number_elements = ((float(num_col)-1)*num_col)/2
	sum = numpy.sum(distance_matrix_array)
	mean_distance = round(sum/number_elements,2)
	max_pairwise_distance = numpy.amax(distance_matrix_array)
	#Transpose the uuper diagonal matrix and convert it into a complete matrix
	#distance_matrix_array_1 = distance_matrix_array.T + distance_matrix_array
	return list_ids,mean_distance,max_pairwise_distance

#Generates the output file (including a header)
def output_generator(output_list,output_path):
	result = []
	header_generator = ["sampled_ID","number_seq","number_seq_without_ambiguities","number_seq_with_ambiguities","aln_length","Variable_sites", "Variable_sites_nuc", "Variable_sites_nuc_parsimony_inf", "Variable_sites_indel","mean_distance","max_pairwise_distance"]
	result.append(header_generator)
	result.append(output_list)

	with open(output_path,"w") as output_path_handle:
		writer = csv.writer(output_path_handle,delimiter ='\t')
		writer.writerows(result)

	return

def wraper_sequence_comparison(input_fasta):
	#output_fasta = alignment(input_fasta,output_fasta)
	input_fasta = convert_to_full_path(input_fasta)
	input_DIR = os.path.dirname(input_fasta)
	aln,number_seq = parse_fasta(input_fasta)
	if number_seq == 1:
		aln_amb_clean,number_seq_without_ambiguities,number_seq_with_ambiguities = check_amb(aln)
		checked_fasta_path = os.path.splitext(input_fasta)[0] + "_checked.fasta"
		checked_output_path = os.path.splitext(checked_fasta_path)[0] + ".out"
		length_aln = wraper_cleaning_alingment(aln_amb_clean,checked_fasta_path)
		sampled_ID = os.path.basename(input_fasta).split(".")[0]
		output_list = [sampled_ID,number_seq,number_seq_without_ambiguities,number_seq_with_ambiguities,length_aln,0,0,0,0,0,0]
		output_generator(output_list,checked_output_path)
		return True
	else:
		aln_amb_clean,number_seq_without_ambiguities,number_seq_with_ambiguities = check_amb(aln)
		checked_fasta_path = os.path.splitext(input_fasta)[0] + "_checked.fasta"
		checked_output_path = os.path.splitext(checked_fasta_path)[0] + ".out"
		length_aln = wraper_cleaning_alingment(aln_amb_clean,checked_fasta_path)
		counter_variable, counter_variable_nuc, counter_variable_nuc_parsimony_inf, counter_variable_indel = seq_comparison(aln_amb_clean)
		path_alignment_fasta_dmat = distmat(checked_fasta_path,input_DIR)
		list_ids,mean_distance,max_pairwise_distance = distmat_to_np_array(path_alignment_fasta_dmat)
	
		#Preparing the output list
		sampled_ID = os.path.basename(input_fasta).split(".")[0]
		output_list = [sampled_ID,number_seq,number_seq_without_ambiguities,number_seq_with_ambiguities,length_aln,counter_variable, counter_variable_nuc, counter_variable_nuc_parsimony_inf, counter_variable_indel,mean_distance,max_pairwise_distance]
		output_generator(output_list,checked_output_path)
		return True

"""
Main
"""

wraper_sequence_comparison(input_fasta)
