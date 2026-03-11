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
parser.add_argument('--input_dir', help='Path of the directory containing subdirectories with the haplotype ID. Every subdirectory must contain only one file with the tree named as concatenate_nuc_gaps.nex_b200.treefile.')
parser.add_argument('--input_hap_info', help='Path to the output with extension .info of the script collapse_haplotypes.py.')
parser.add_argument('--input_reference_clade', help='Path of the file containing the information reference tree. It must be a file with two columns separated by tab. The fisrt columns must contain the clade name (arbitrary one word ID) and the second the list of the samples corresponding to the clade, separated by comma.')
parser.add_argument('--output', help='A tab delimited file to store the results. It contains the following fields: samplesID, haplotypeID and the information for the clades.')


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()

if args.input_dir:
	input_dir = args.input_dir
	if input_dir[-1] != "/":
		input_dir = input_dir + "/"

if args.output:
	output = args.output

if args.input_hap_info:
	input_hap_info = args.input_hap_info

if args.input_reference_clade:
	input_reference_clade = args.input_reference_clade

'''
FUNCTIONS
'''

def convert_to_full_path(path):
	#Get working directory
	absolute_path = os.path.abspath(path)
	return absolute_path

#Parse the file contaning the haplotype information
def parse_ref_hap_info(ref_hap_info):
	ref_hap_info_dict = {}
	with open(ref_hap_info) as file:
		for line in file:
			info = line.rstrip()
			info_list = info.split("\t")
			ref_hap_info_dict[info_list[0]] = info_list[2].split(",")

	return ref_hap_info_dict

#Parse the file contaning the information of the reference phylogeny
def file_to_list(input_file):
	clades_dict = {}
	with open(input_file) as file:
		all_ids_ref_clades = []
		for line in file:
			info = line.rstrip()
			info_list = list(info.split("\t"))
			clades_dict[info_list[0]] = list(info_list[1].split(","))
			for name in info_list[1].split(","):
				all_ids_ref_clades.append(name)
	all_ids_ref_clades.sort()
	return clades_dict,all_ids_ref_clades

#Obtin the ids of a tree
def get_ids_from_trees(tree_file):
	ref_tree = Tree(tree_file, format=0)
	all_ids_current_tree = []
	for node in ref_tree.traverse("postorder"):
		if node.is_leaf():
			all_ids_current_tree.append(node.name)
	all_ids_current_tree.sort()
	return all_ids_current_tree

#Main function to produce the results
def comparing_ids(all_ids_ref_clades,all_ids_current_tree,current_sample,current_general_id,clades_dict,current_tree):
	clades_dict_current = copy.deepcopy(clades_dict)
	all_ids_plus_sample = copy.deepcopy(all_ids_ref_clades)
	all_ids = copy.deepcopy(all_ids_ref_clades)
	all_ids_plus_sample.append(current_sample)
	all_ids_plus_sample.sort()

	#Full phylogeny. All reference sequences + current sequence
	if all_ids_plus_sample == all_ids_current_tree:
		for clade in clades_dict_current.keys():
			clades_dict_current[clade].append(current_sample)
		
		#Result of the new intraspecific haplotype
		current_result_new_intraspecific_haplotype,clade_result_new_intraspecific_haplotype = monophyly_test(clades_dict_current,current_tree)		
		
		current_result_most_coverage,clade_result_most_coverage = get_haplotype_info(ref_hap_info_dict,current_general_id,clades_dict)

		return clades_dict_current,current_result_new_intraspecific_haplotype,current_result_most_coverage,clade_result_new_intraspecific_haplotype,clade_result_most_coverage

	#Edited phylogeny. ID replacement. All reference sequences without interspecies haplotype of the most covered intresp haplotype + current sequence
	elif len(all_ids) == len(all_ids_current_tree):
		ids_only_current_tree = []
		ids_only_ref_tree = []
		for id in all_ids:
			if id not in all_ids_current_tree:
				ids_only_ref_tree.append(id)
		for id in all_ids_current_tree:
			if id not in all_ids:
				ids_only_current_tree.append(id)
		
		for clade in clades_dict_current.keys():
			for id in ids_only_ref_tree:
				if id in clades_dict_current[clade]:
					current_clade = clades_dict_current[clade]
					current_clade.remove(id)
					clades_dict_current[clade] = current_clade
			
			current_clade = clades_dict_current[clade]
			current_clade.append(current_sample)
			clades_dict_current[clade] = current_clade
		
		#Result of the new intraspecific haplotype
		current_result_new_intraspecific_haplotype,clade_result_new_intraspecific_haplotype = monophyly_test(clades_dict_current,current_tree)			
		current_result_most_coverage,clade_result_most_coverage = get_haplotype_info(ref_hap_info_dict,current_general_id,clades_dict)

		return clades_dict_current,current_result_new_intraspecific_haplotype,current_result_most_coverage,clade_result_new_intraspecific_haplotype,clade_result_most_coverage
		#clades_dict_current,ids_only_ref_tree,ids_only_current_tree
	else:
		print("something goes wrong in comparing_ids function 1")

	

	
#Test for monophyly using a dictionary of node names as key and leave names as element.
def monophyly_test(selected_leaves,tree):
	t = Tree(tree)
	t.set_outgroup( t&"Agrandis_8.129" )
	current_result = []
	#selected_leaves_copy = dict(selected_leaves)
	for test_group in selected_leaves:
		current_result.append(test_group)
		if t.check_monophyly(values=selected_leaves[test_group], target_attr="name")[0] == True:
			current_result.append("yes")
			clade_result = test_group
		else:
			current_result.append("no")
	return current_result,clade_result

def get_haplotype_info(ref_hap_info_dict,current_general_id,clades_dict):
	id_ref_tree = []
	current_result = []
	for ids_hap in ref_hap_info_dict.keys():
		for id in ref_hap_info_dict[ids_hap]:
			if id.startswith(current_general_id):
				id_ref_tree.append(ref_hap_info_dict[ids_hap][0])
	#print(clades_dict_current)
	for clade in clades_dict.keys():
		current_result.append(clade)
		tmp = []
		for id in clades_dict[clade]:
			if id_ref_tree[0] == id:
				tmp.append("yes")
				break
		if tmp == ['yes']:
			current_result.append('yes')
			clade_result = clade
		else:
			current_result.append('no')

	return current_result,clade_result

def save_results(final_results,output):
	with open(output,"w") as output_file:
		for current_result in final_results:
			current_result_new_intraspecific_haplotype = "". join(current_result[2])
			current_result_most_coverage = current_result[3]
			output_file.write(current_result[0] + "\t" + current_result[1] + "\t" + "\t". join(current_result[2]) + "\t" + "\t". join(current_result[3]) + "\t" + current_result[4] + "\t" + current_result[5]   + "\n")
	return

'''
Main
'''

final_results = []
input_dir = convert_to_full_path(input_dir)
ref_hap_info = convert_to_full_path(input_hap_info)
input_file = convert_to_full_path(input_reference_clade)
ref_hap_info_dict = parse_ref_hap_info(ref_hap_info)
clades_dict,all_ids_ref_clades = file_to_list(input_file)

for current_sample in os.listdir(input_dir):
	current_general_id = ("_").join(current_sample.split("_")[:-1])
	sampleDIR = os.path.join(input_dir,current_sample)
	if os.path.exists(os.path.join(sampleDIR,"concatenate_nuc_gaps.nex_b200.treefile")):
		current_tree = os.path.join(sampleDIR,"concatenate_nuc_gaps.nex_b200.treefile")
		all_ids_current_tree = get_ids_from_trees(current_tree)
		current_result = []
		clades_dict_current,current_result_new_intraspecific_haplotype,current_result_most_coverage,clade_result_new_intraspecific_haplotype,clade_result_most_coverage = comparing_ids(all_ids_ref_clades,all_ids_current_tree,current_sample,current_general_id,clades_dict,current_tree)
		current_result.append(current_general_id)
		current_result.append(current_sample)
		current_result.append(current_result_new_intraspecific_haplotype)
		current_result.append(current_result_most_coverage)
		current_result.append(clade_result_new_intraspecific_haplotype)
		current_result.append(clade_result_most_coverage)
		final_results.append(current_result)
	elif os.path.exists(os.path.join(sampleDIR,"concatenate_nuc_gaps.nex_b200_edited.treefile")):
		current_tree = os.path.join(sampleDIR,"concatenate_nuc_gaps.nex_b200_edited.treefile")
		all_ids_current_tree = get_ids_from_trees(current_tree)
		current_result = []
		clades_dict_current,current_result_new_intraspecific_haplotype,current_result_most_coverage,clade_result_new_intraspecific_haplotype,clade_result_most_coverage = comparing_ids(all_ids_ref_clades,all_ids_current_tree,current_sample,current_general_id,clades_dict,current_tree)
		current_result.append(current_general_id)
		current_result.append(current_sample)
		current_result.append(current_result_new_intraspecific_haplotype)
		current_result.append(current_result_most_coverage)
		current_result.append(clade_result_new_intraspecific_haplotype)
		current_result.append(clade_result_most_coverage)
		final_results.append(current_result)
	else:
		print("There is no a phylogeny file in " + sampleID)


save_results(final_results,output)

