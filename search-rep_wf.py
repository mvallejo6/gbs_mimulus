from Bio import SeqIO
from Bio import Seq
import re 

#python program location:

datapath = "../genome_data/"
#filename = "mimseq.fasta"
#original_file = "mimseq.fasta"
original_file = "Mimulus_guttatus_luteus_combined.fasta"
corrected_file = "m_gut_lut_corrected.fasta"
scaffolds_file = "scaffolds.txt"

with open(datapath + original_file) as original, open(datapath + corrected_file, 'w') as corrected, open(datapath + scaffolds_file, 'w') as scaffolds:
    records = SeqIO.parse(datapath + original_file, 'fasta')

    for record in records:
	    print("Sequence Name: %s, length: %i" % (record.id, len(record)))
            record.id = record.id.replace('chr','gut_')
            record.id = record.id.replace('lcl|scaffold_','lut_')
            record.description = ""
            print("Sequence Name: %s, length: %i" % (record.id, len(record)))
            SeqIO.write(record, corrected, 'fasta')
            scaffolds.write(record.id + "," + str(len(record)) + '\n')
            
