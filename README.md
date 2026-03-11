# Intragenomic variation

This repository contains the scripts used in the paper:

> **Single-molecule amplicon sequencing reveals phylogenetic consistency across intraindividual ITS1 copies of the fruit fly, *Anastrepha fraterculus* complex (Diptera: Tephritidae)**  
> *[Add full citation here: authors, year, journal, volume, pages, DOI]*

---

## 1. alignment_variation.py
This script calculates the following variation parameters from a multifasta file: 

1. Number of sequences
2. Number of sequences without ambiguities
3. Number of sequences with ambiguities
4. Alignment length 
5. Number of variable sites including nucleotides and indels
6. Number of variable sites including only nucleotides
7. Number of parsimony informative sites  
8. Number of variable sites including only indels
9. Mean of pairwise distance including gaps (gapweight=0.5) calculated using distmat tool (https://www.bioinformatics.nl/cgi-bin/emboss/distmat)
10. Maximum pairwise distance including gaps (gapweight=0.5)  calculated using distmat tool (https://www.bioinformatics.nl/cgi-bin/emboss/distmat)


## 1.1 Requirements

The script requires:
- mafft (https://mafft.cbrc.jp/alignment/software/)
- distmat (https://emboss.sourceforge.net/apps/cvs/emboss/apps/distmat.html)
- Python 3 and the following modules:
  - Biopython
  - numpy

## 1.2 Usage
```bash
python alignment_variation.py  --input_fasta input.fasta
```

---


## 2. collapse_haplotypes.py
This script collapses haplotypes. 
Important:
This script does not support "N" or ambiguous bases in the input sequences.
This script separates sequences on haplotypes based on nucleotides and indels.


## 2.1 Requirements

The script requires:
- Python 3 and the following modules:
  - Biopython

## 2.2 Usage
```bash
python  collapse_haplotypes.py  --input_fasta  input.fasta --output_info  output.inf --output_fasta output.fasta
```

---


## 3. retrieve_fasta_coordinates.py
Retrieve fasta sequences based on a list of ids and coordinates from a multifasta file (database).
Important:
If the first coordinate is greater than the second coordinate in the input_list file, the program will generate a reverse complement sequence.


## 3.1 Requirements

The script requires:
- Python 3 and the following modules:
  - Biopython

## 3.2 Usage
```bash
python retrieve_fasta_dict.py --input_list list_IDs_coordinates  --input_fasta input.fasta --output_fasta output.fasta --mode fast
```

---


## 4. analyze_trees.py
Retrieve fasta sequences based on a list of ids and coordinates from a multifasta file (database).
Important:
If the first coordinate is greater than the second coordinate in the input_list file, the program will generate a reverse complement sequence.


## 4.1 Requirements

The script requires:
- Python 3 and the following modules:
  - ETE3

## 4.2 Usage
```bash
python analyze_trees.py --input_dir haplotypes --input_hap_info input_hap.info --input_reference_clade phylo_groups --output output.log
```

---

## Reference 

If you use these scripts, please cite the original paper:

> **Single-molecule amplicon sequencing reveals phylogenetic consistency across intraindividual ITS1 copies of the fruit fly, *Anastrepha fraterculus* complex (Diptera: Tephritidae)**  
> *[Add full citation here]*

This script is in the public domain in the United States per 17 U.S.C. § 105
