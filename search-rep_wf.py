from Bio import SeqIO
from Bio import Seq
import re 

#This script changes the names of the scaffolds in the merged guttayus-luteus genome

datapath = "../genome_data/"
#filename = "mimseq.fasta"
#original_file = "mimseq.fasta"

#File with the original data in fasta format
original_file = "Mimulus_guttatus_luteus_combined.fasta"
#Name of the file where new results will be stored
corrected_file = "m_gut_lut_corrected.fasta"
# A log file that will be created to list all the scaffold names and length of their sequence
scaffolds_file = "scaffolds.txt"

#The "with" command allows open and automatically closing files. It also allows you to 
#open and define file names in a single command
with open(datapath + original_file) as original, open(datapath + corrected_file, 'w') as corrected, open(datapath + scaffolds_file, 'w') as scaffolds:
    records = SeqIO.parse(datapath + original_file, 'fasta')

    for record in records:
	    print("Sequence Name: %s, length: %i" % (record.id, len(record))) #Prints to screen
            #M. guttatus records
            record.id = record.id.replace('chr','gut_') #Replaces the 'chr' string in the record.id field with 'gut_'
            #M. luteus records
            record.id = record.id.replace('lcl|scaffold_','lut_')
            #Delete the record description (if present)
            record.description = ""
            print("Sequence Name: %s, length: %i" % (record.id, len(record)))
            #Write the new record to a new file
            SeqIO.write(record, corrected, 'fasta')
            #Write the log file with scaffold names and sequence length
            scaffolds.write(record.id + "," + str(len(record)) + '\n')
            
