import pandas as pd 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# PflA aln length 924
# PflB aln length 987

pfla_d = pd.read_excel('../../1_bidirectional-best-hits/results-CjPflA-names.xlsx', sheet_name=['BBH-eval-filter', 'psiblast-eval-filter'])
pflb_d = pd.read_excel('../../1_bidirectional-best-hits/results-CjPflB-names.xlsx', sheet_name=['BBH-eval-filter', 'psiblast-eval-filter'])

# make a list of tuples, first element is accver of hit protein, second is gcf of that organism
def acc_gcf(dic):
	bbh = dic['BBH-eval-filter']
	psi = dic['psiblast-eval-filter']
	pfl_bbh = [(bbh.iloc[i,1], bbh.iloc[i,9]) for i in range(len(bbh))]
	pfl_psi = [(psi.iloc[i,1], psi.iloc[i,9]) for i in range(len(psi))]
	pfl = pfl_bbh+pfl_psi
	return pfl

pfla = acc_gcf(pfla_d)
pflb = acc_gcf(pflb_d)

with open('PflAB-paired-alignment.fasta', 'a+') as outfile:
	for record in SeqIO.parse('PflA-protein-nr-ids.fasta', 'fasta'):
		acc_a = record.id
		desc = record.description.split('[')[1][:-1]
		for ela in pfla: # go through the pflA hits and find the one that matches the accession
			if acc_a in ela[0]:
				gcf = ela[1] # get the gcf of the match 
				for elb in pflb:
					if gcf == elb[1]: # find that gcf in pflb hits
						acc_b = elb[0]
						for rec in SeqIO.parse('PflB-protein-nr-ids.fasta', 'fasta'):
							if acc_b in rec.id:
								cat_seq = record.seq + rec.seq
								r = SeqRecord(cat_seq, id=gcf, description=desc)
								SeqIO.write(r, outfile, 'fasta')
					else: 
						print('not in pflB:', gcf)
