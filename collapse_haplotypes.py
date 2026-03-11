#!/usr/bin/env python

"""
This script collapses haplotypes. 
Important:
This script does not support "N" or ambiguous bases.
This script separates sequences on haplotypes based on nucletides and indels.
"""

from Bio import SeqIO
import itertools
import argparse,os,sys

class MyParser(argparse.ArgumentParser):
        def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_fasta', help='Path to the alignment in fasta format.')
parser.add_argument('--output_info', help='Tab delimited file informing id of one of the sequences of the haplotype, frequency and ids of all sequences of the haplotype.')
parser.add_argument('--output_fasta', help='Path to the collapsed alignment.')

if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

args = parser.parse_args()


if args.input_fasta:
        input_fasta = args.input_fasta

if args.output_fasta:
        output_fasta = args.output_fasta
        
if args.output_info:
        output_info = args.output_info

"""
Main
"""

records = list(SeqIO.parse(input_fasta, "fasta"))

dict_unique = {}
list_record_ids = []
for record in records:
        list_record_ids.append(record.id)

if len(records) > 1:
        total_number_combination = ((len(records)-1)*len(records))/2
        counter = 0
        for record1,record2 in itertools.combinations(records, 2):
                if record1.id not in dict_unique.keys():
                        dict_unique[record1.id]=[record1.id]
                if record2.id != None:
                        if str(record1.seq) == str(record2.seq):
                                dict_unique[record1.id].append(record2.id)
                                index = list_record_ids.index(record2.id)
                                records[index].id = None
                counter = counter + 1
                if counter == total_number_combination:
                        if record2.id != None:
                                dict_unique[record2.id]=[record2.id]

        #Clean the dictionary
        if None in dict_unique.keys():
                del dict_unique[None]

elif len(records) == 1:
        dict_unique[records[0].id]=[records[0].id]
        

#Retrieve fasta and info
record_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
with open(output_fasta, "w") as output_handle:
        with open(output_info, "w") as output_info_handle:
                for ID_seq_record in dict_unique.keys():
                        if ID_seq_record in record_dict:
                                #Save the fasta
                                output_handle.write(record_dict[ID_seq_record].format("fasta"))
                                #Save the haplotype info
                                output_info_handle.write(ID_seq_record + "\t" + str(len(dict_unique[ID_seq_record])) + "\t" + ",".join(dict_unique[ID_seq_record]) + "\n")
                        else:
                                print(ID_seq_record + " is not in " + input_fasta)

