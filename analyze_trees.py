#!/usr/bin/env python

import argparse,os,sys
import copy
from ete3 import Tree


class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', help='Path of the directory containing subdirectories with the haplotype ID. Every subdirectory must contain two files sampleID_list and the tree named as concatenate_nuc_gaps.nex_b200.treefile.')
parser.add_argument('--input_reference_clade', help='File with reference haplotypes and clade name in a tab delimited file.')
parser.add_argument('--outgroup', help='ID of the outgroup as is in the tree')
parser.add_argument('--output', help='A tab delimited file to store the results. It contains the following fields: HaplotypeID, clade of the dominant haplotype, and clade of the alternate haplotype.')


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()

if args.input_dir:
	input_dir = args.input_dir
	if input_dir[-1] != "/":
		input_dir = input_dir + "/"

if args.input_reference_clade:
	input_reference_clade = args.input_reference_clade

if args.outgroup:
	outgroup = args.outgroup

if args.output:
	output = args.output


'''
FUNCTIONS
'''

def convert_to_full_path(path):
	#Get working directory
	absolute_path = os.path.abspath(path)
	return absolute_path

def process_input_reference_clade(input_reference_clade):
	reference_dict = {}
	with open(input_reference_clade, 'r') as f:
		for line in f:
			current_ref_hap_clade = line.rstrip().split("\t")
			reference_dict[current_ref_hap_clade[0]] = current_ref_hap_clade[1]
	return reference_dict

def current_reference_list_to_dictionary_clade(reference_dict,current_reference_list,sampleID,input_reference_clade):
	reference_clade_association_dict = {}
	
	for current_ref_hap in current_reference_list:
		if current_ref_hap in reference_dict:
			clade = reference_dict[current_ref_hap]
			if clade in reference_clade_association_dict.keys():
				reference_clade_association_dict[clade].append(current_ref_hap)
			else:
				reference_clade_association_dict[clade] = []
				reference_clade_association_dict[clade].append(current_ref_hap)

		else:
			print("Check " + current_ref_hap + ". This ID is not in " + input_reference_clade)

	return reference_clade_association_dict

def get_current_reference_list(input_list,sampleID):
	current_reference_list = []
	with open(input_list, 'r') as f:
		for line in f:
			if line.rstrip() != sampleID:
				current_reference_list.append(line.rstrip())

	return current_reference_list

def get_expected_clade(current_reference_list,reference_dict):
	reference_list = []
	for ref_hap_id in reference_dict.keys():
		reference_list.append(ref_hap_id)
	reference_list_set = set(reference_list)
	current_reference_list_set = set(current_reference_list)
	expected_ID = list(reference_list_set - current_reference_list_set)[0]
	expected_clade = reference_dict[expected_ID]
	return expected_clade

#Test for monophyly using a dictionary of node names as key and leave names as element.
def monophyly_test(reference_clade_association_dict,tree,outgroup,sampleID,expected_clade):
	t = Tree(tree)
	t.set_outgroup( t&outgroup )
	#Just in case something gets wrong
	current_result = [sampleID,expected_clade,"NA"]

	#selected_leaves_copy = dict(selected_leaves)
	for test_group in reference_clade_association_dict.keys():

		test_group_list = reference_clade_association_dict[test_group]
		test_group_list.append(sampleID)
		if t.check_monophyly(values=test_group_list, target_attr="name")[0] == True:
			current_result = [sampleID,expected_clade,test_group]
			break

	return current_result

def save_results_file(final_result,output):
	header = ["SampleID","Dominant_clade","Alternate_clade"]
	with open(output,"w") as out_f:
		out_f.write(("\t").join(header) + "\n")
		for current_result in final_result:
			out_f.write(("\t").join(current_result) + "\n")
	return


'''
Main
'''

input_dir = convert_to_full_path(input_dir)
input_reference_clade = convert_to_full_path(input_reference_clade)
reference_dict = process_input_reference_clade(input_reference_clade)

final_result = []

for subdir in os.listdir(input_dir):
	sampleID = subdir
	list_file =  os.path.join(input_dir,subdir,sampleID + "_list")
	tree_file =  os.path.join(input_dir,subdir,"concatenate_nuc_gaps.nex_b200.treefile")
	#Get only reference
	current_reference_list = get_current_reference_list(list_file,sampleID)
	#Get expected clade (dominant clade)
	expected_clade = get_expected_clade(current_reference_list,reference_dict)
	#Get dictionary clade and leaves in the clades to test monophyly
	reference_clade_association_dict = current_reference_list_to_dictionary_clade(reference_dict,current_reference_list,sampleID,input_reference_clade)
	#Test monophyly and get results
	current_result = monophyly_test(reference_clade_association_dict,tree_file,outgroup,sampleID,expected_clade)
	#Save the results into a list
	final_result.append(current_result)

#Save the result into a file
save_results_file(final_result,output)
