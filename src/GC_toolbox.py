#!/usr/bin/python3

##########################################################################################################################
#PROGRAM NAME   = GC_toolbox.py
#DESCRIPTION    = Selectively parses whole genome or CDS fasta/fastq and produces overview table
#                 + plots for mean GC content across sequences and, in case of CDS, codon positions
#INPUT FILES    = Any fastq or DNA fasta file
#OUTPUT FILES   = Converts fastq to a fasta file. Produces both .tsv and .pdf overview table + two pdf files for GC plots
#RUN BY         = python3 GC_toolbox.py input_file [whole_genome]/[CDS]
#AUTHOR 	    = Federico Musso
#VERSION 1.0	14.03.2020
##########################################################################################################################

import sys

#define plastic function to count GCs which will be called in the following
#function
def gc_counter(seq):
    gc = seq.count('G') + seq.count('C')
    stdnuc = seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T')
    GC_percent = round(100* gc/stdnuc,1)
    return(GC_percent)
#define function for GC calculation that will operate on a sequence variable.
def gc_calculations(sequences):
    gc = gc_counter(sequences)          #we call the previous function on this one
    first_pos = sequences[::3]          #start from beginning of the string and returns every 3rd pos
    gc_first = gc_counter(first_pos)    #then we apply the count function
    second_pos = sequences[1::3]        #starts from second position, same behaviour
    gc_second = gc_counter(second_pos)
    third_pos = sequences[2::3]
    gc_third = gc_counter(third_pos)
    return('{0}\n{1}\n{2}\n{3}\n'.format(gc, gc_first, gc_second, gc_third)) #return all four GC calculations on separate lines
#define function for reverse complement
def rev_complementer(sequence):
    rev = ""
    for i in sequence:
        rev = i + rev
    revcomp = rev.translate(str.maketrans('ACGT','TGCA')) #see compendium
    return revcomp

# # to test the function with a toy sequence:
# funct_out= gc_calculations('AATGANACGATGATCCA')
# print(funct_out)  #we have to print the output since we used return and not directly print in the function
# # we feed the sequence to the function -> the sequence is stored into the "sequences"
# # variable and we can then perform various operations on the latter



# we open two files. Outfile will be used for the production of a Fasta file

#########################################################################################################################
#IDEA: parse and operate on FASTQ or FASTA file when the input is a whole genome, parse FASTA file when dealing with CDS
#########################################################################################################################

#start reading file
with open(sys.argv[1], 'r') as fin:
    sequencetype = sys.argv[2]          #introduce whole_genome or CDS. Second argument will make code execution selective.
    fin = fin.readlines()               #read file line by line

    #initiate empty lists which will contain sequence lengths
    lengths_list = []
    GcTot_list = []
    GcFirst_list = []
    GcSecond_list = []
    GcThird_list = []
    revcomp_list = []
    #initiate empty dictionary for sequences and IDs
    sequencedict = {}

    #################
    # WHOLE GENOME  #
    #################
    #conditional execution for whole_genome input file, parses a FASTQ or a FASTA file
    if sequencetype == "whole_genome":
        for line in fin:                #here we save the first line only and use it later to evaluate whether
            firstline = line            #we are dealing with a fasta or fastq file
            break

        #FASTA parser
        if firstline.startswith(">"):   #if we had a fasta file, start treating it as such and create dictionary
            for line in fin:
                if line.startswith(">"):
                    seqmemory = ''
                    current_ID = line
                    sequencedict.update({current_ID:''})
                else:
                    seqmemory = seqmemory + line
                    seqmemory = seqmemory.rstrip("\n")
                    sequencedict[current_ID] = seqmemory

        #FASTQ parser
        elif firstline.startswith("@"):
            with open('./results/out.fasta', 'w') as fout:                 #output file is opened only if dealing with FASTQ
                for count,line in enumerate(fin):
                    if count % 4 == 0:                              #isolate ID lines
                        current_ID = line.rstrip()
                        current_ID = current_ID.split(" ")[0][1:]   #isolate the proper ID characters
                        sequencedict.update({current_ID:''})        #set ID as key
                    if count % 4 == 1:                              #isolate sequences
                        seq = line.rstrip()
                        seq = seq.upper()
                        sequencedict[current_ID] = seq              #put sequence as a content for that ID
                        print(">"+current_ID+"\n"+seq, file=fout)       #conversion to FASTA file
        else:
            sys.exit("\nERROR!! Your file is neither a FASTA nor a FASTQ file!\n")
        #now we have a dictionary with IDs as keys and sequences as content. we loop through it!
        for key in sequencedict:
            seq = sequencedict[key]
            #whole genome: perform mean GC calculations and create FASTA file
            totGC = float(gc_calculations(seq)[0:4])
            GcTot_list.append(totGC)

    ###################
    # CDS FASTA input #
    ###################
    #if we specified that we are dealing with a CDS, it will parse a FASTA file
    if sequencetype == "CDS":
        seqmemory = ''
        oksymbols = ('A','C','T','G','N')  #add tolerated nucleotide characters. Only works with tuples of strings
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                seqmemory = ''
                current_ID = line
                sequencedict.update({current_ID:''})
            elif not line.startswith(oksymbols):  #check for okay symbols
                sys.exit("\n\nERROR! Input file doesn't match specified sequence type: might be a FASTQ, or there's some weird sequence!\n\n")
            else:
                seqmemory = seqmemory + line
                seqmemory = seqmemory.rstrip("\n")
                sequencedict[current_ID] = seqmemory

        #we use the same sequence dictionary to perform the calculation.
        for key in sequencedict:
            seq = sequencedict[key]
            totGC = gc_calculations(seq)[0:4]              #this enables us to retrieve only some
            totGC = float(totGC)                           #characters (=the first GC value of the four calculations)
            GcTot_list.append(totGC)                       #from the function output. These are stored in a list

            firstGC = float(gc_calculations(seq)[5:9])     #isolate GC1 values. Smoother syntax for floating
            GcFirst_list.append(firstGC)

            secondGC = float(gc_calculations(seq)[10:14])  #isolate GC2 values
            GcSecond_list.append(secondGC)

            thirdGC = float(gc_calculations(seq)[15:19])   #isolate GC3 values
            GcThird_list.append(thirdGC)

        #no conversion to FASTA since we already gave a FASTA file as input

    ######################
    # General processing #
    ######################

    # General operations on the same dictionary, to produce a summary table
    for key in sequencedict:
        seq = sequencedict[key]
        length = len(seq)                      #produce list with sequence length
        lengths_list.append(length)
        revcomp = rev_complementer(seq)        #Produce REVERSE COMPLEMENT list, for both files
        revcomp_list.append(revcomp)

    #MAX,MIN,MEAN LENGTH CALCULATION AND OVERALL GC CONTENT
    #iterate through list containing sequence lengths and extract mean length
    maxlength = max(lengths_list)
    minlength = min(lengths_list)
    tot = len(lengths_list)
    sum = 0
    for i in lengths_list:
        sum += i
    mean_length = round(sum/tot,2)

    #iterate through list containing GC percentages and extract mean GC%.
    #NOT A LIST, just a single value to display on a header table
    totGC = len(GcTot_list)
    GCsum = 0
    for i in GcTot_list:
        GCsum += i
    mean_GC = round(GCsum/totGC,2)

#PRODUCE OVERVIEW TABLE WITH SEQUENCE INFO
#first as a .tsv file
with open('./results/summary_table.tsv', 'w') as fout:
    print("Input file\t"+sys.argv[1]+"\n"+"Data type\t"+sequencetype+"\n"+ \
    "Tot sequences\t"+str(tot)+"\n"+"Mean GC%\t"+str(mean_GC)+"\n"+"Average length\t"+ \
    str(mean_length)+"\n"+"Longest seq\t"+str(maxlength)+"\n"+"Shortest seq\t"+str(minlength)+"\n",file=fout)

#then process with pandas
import matplotlib.pyplot as plt
import pandas
from pandas.plotting import table

df = pandas.read_csv("./results/summary_table.tsv", sep="\t", header=None) #read into dataframe
#df.reset_index(drop=True, inplace=True)
df.index = ['']*7 #remove horizontal indexes
df.columns = ['','Observed parameters']
ax = plt.subplot(721, frame_on=False) # no visible frame
ax.xaxis.set_visible(False)  # hide the x axis
ax.yaxis.set_visible(False)  # hide the y axis

table(ax,df)
plt.savefig('./results/summary_table.pdf')
plt.close()

# PLOTTING OF THE OVERALL PER SEQUENCE GC CONTENT
# done with matplotlib, both for whole genomes and CDS
plt.hist(GcTot_list,edgecolor= 'black', bins = 30, rwidth = 0.95, label = 'Mean GC content')
plt.xlabel('GC percentage')
plt.ylabel('Observed frequency')
plt.title('Mean GC content across sequences')
plt.legend()
plt.savefig("./results/TotalGC.pdf")
plt.close()
# PLOTTING OF THE GC CONTENT IN FIRST, SECOND AND THIRD POSITION
# only for CDS datatype
if sequencetype == "CDS":
    import seaborn as sns

    # plt.hist([GcFirst_list,GcSecond_list,GcThird_list], color= ['red','green','cyan'],
    # edgecolor= 'black', bins= 20, rwidth=0.90, label=['GC1','GC2','GC3'])
    sns.distplot(GcFirst_list, hist= False, kde=True,               #IMPROVE: avoid code repetition
    color = 'orange', kde_kws = {'shade': True, 'linewidth': 2},
    label = 'GC1')
    sns.distplot(GcSecond_list, hist= False, kde=True,
    color = 'cyan', kde_kws = {'shade': True, 'linewidth': 2},
    label = 'GC2')
    sns.distplot(GcThird_list, hist= False, kde=True,
    color = 'green', kde_kws = {'shade': True, 'linewidth': 2},
    label = 'GC3')
    plt.xlabel('GC percentage')
    plt.ylabel('Observed frequency')
    plt.title('GC content in first, second and third codon position')
    plt.legend(prop={'size': 8}, title = 'Codon positions')
    plt.savefig("./results/CodonGC.pdf")
    plt.close()
