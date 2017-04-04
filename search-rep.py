from Bio import SeqIO
from Bio import Seq
import re 

#python program location: cd ~/projects/mim-uv/analyses
#datapath = "../../mluteus_genome/"
filename = "mimseq.fasta"

###To read FASTA files use:
#records = SeqIO.parse(datapath + filename, "fasta")
records = SeqIO.parse(filename, "fasta")
# records behaves like a list

#for target in records:
#	print("Sequence ID: %s, length: %i" % (target.id, len(target)))  #Usually one record only

for target in records:
	print("Sequence Name: %s, length: %i" % (target.id, len(target)))
	#print("New record id: %s" target.id.replace("gb", "NEW"))
	print(target.id)
	print(target.name)
	temp = str(target.name)
	print(temp)
	temp = temp.replace("|gb|", "NEW")
	#temp.rep
	print("New name: %s" % str(temp))
