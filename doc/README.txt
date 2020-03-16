GC_toolbox.py
V1.0		Federico Musso 		14.03.2020

This program is a simple python script for extracting 
sequence information from FASTQ and FASTA files. It 
produces a summary table with some key information on 
the sequences contained in the input file, including: 
total number of sequences in the input file, mean GC%, 
average sequence length, longest and shortest sequence 
in the input file. Additionally, it plots per sequence 
GC content and returns it as a pdf file. If the input 
data represents a CDS, GC content for each codon 
position will be plotted, too.

#INSTALLATION
Clone from Github or use the package manager

#USAGE
The script is called in the following way:

python3 GC_toolbox.py input_file [whole_genome]/[CDS]

note: relative pathways in the script are set for 
storage of results in a "results" folder. Paths should 
be changed in case this is not preferred. 

#CONTRIBUTING
Open for any contribution


