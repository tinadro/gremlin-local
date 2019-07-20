from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

with open('PflA.afa', 'r') as a, open('PflB.afa', 'r') as b, open('PflAB-joined.afa', 'a+') as out: # second sys arg names output file
	seq_desc = list()
	for seq in SeqIO.parse(a, 'fasta'): # first sys arg names the file you wanna cleanup 
		print(seq.id)
		#for seq2 in SeqIO.parse(b, 'fasta'):
			#print(seq.description)
			#if seq.description in seq2.description:
			#	sequence = seq.seq + seq2.seq
			#	record = SeqRecord(sequence, id=seq.id, name=seq.name, description=seq.description)
			#	print(seq.id, seq.name, seq.description)
			#	SeqIO.write(record, out, 'fasta')
				
#		if seq.description 
#		if seq.id not in seq_ids: #if the seqid has not appeared in the file yet
#			seq_ids.append(seq.id)
#			SeqIO.write(seq, outfile,'fasta') # add it to the output file. otherwise do nothing. 

