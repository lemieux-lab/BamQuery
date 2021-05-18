
# Translation tables
# Ref: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
translTable = dict()

# Standard Code (NCBI transl_table=1)
translTable['default'] = {
'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',

'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M',
'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T',
'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K',
'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R',

'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G',

'!GA' : 'U'

}
codonTable = translTable['default']


# The Vertebrate Mitochondrial Code (NCBI transl_table=2)
translTable['mt'] = {
'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : 'W',
'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',

'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'M', 'ATG' : 'M',
'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T',
'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K',
'AGT' : 'S', 'AGC' : 'S', 'AGA' : '*', 'AGG' : '*',

'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G'
}

def complementTab(seq=[]):
	"""returns a list of complementary sequence without inversing it"""
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
				  'W': 'W', 'S': 'S', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', 'a': 't',
				  'c': 'g', 'g': 'c', 't': 'a', 'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm', 'w': 'w',
				  's': 's', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n'}
	seq_tmp = []
	for bps in seq:
		if len(bps) == 0:
			#Need manage '' for deletion
			seq_tmp.append('') 
		elif len(bps) == 1:
			seq_tmp.append(complement[bps])
		else:
			#Need manage 'ACT' for insertion
			#The insertion need to be reverse complement (like seq)
			seq_tmp.append(reverseComplement(bps))
			
	#Doesn't work in the second for when bps==''
	#seq = [complement[bp] if bp != '' else '' for bps in seq for bp in bps]

	return seq_tmp
	
def reverseComplementTab(seq):
	'''
	Complements a DNA sequence, returning the reverse complement in a list to manage INDEL.
	'''
	return complementTab(seq[::-1])

def reverseComplement(seq):
	'''
	Complements a DNA sequence, returning the reverse complement.
	'''
	return complement(seq)[::-1]

def complement(seq) :
	"""returns the complementary sequence without inversing it"""
	tb = str.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
							  "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")

	#just to be sure that seq isn't unicode
	return str(seq).translate(tb)

def translateDNA(sequence, frame = 'f1', translTable_id='default') :
	"""Translates DNA code, frame : fwd1, fwd2, fwd3, rev1, rev2, rev3"""

	protein = ""

	if frame == 'f1' :
		dna = sequence
	elif frame == 'f2':
		dna = sequence[1:]
	elif frame == 'f3' :
		dna = sequence[2:]
	elif frame == 'r1' :
		dna = reverseComplement(sequence)
	elif frame == 'r2' :
		dna = reverseComplement(sequence)
		dna = dna[1:]
	elif frame == 'r3' :
		dna = reverseComplement(sequence)
		dna = dna[2:]
	else :
		raise ValueError('unknown reading frame: %s, should be one of the following: fwd1, fwd2, fwd3, rev1, rev2, rev3' % frame)

	for i in range(0, len(dna),  3) :
		codon = dna[i:i+3]

		# Check if variant messed with selenocysteine codon
		if '!' in codon and codon != '!GA':
			codon = codon.replace('!', 'T')

		if (len(codon) == 3) :
			try :
				# MC
				protein += translTable[translTable_id][codon]
			except KeyError :
				print ('sequence ',sequence)
				combinaisons = polymorphicCodonCombinaisons(list(codon))
				translations = set()
				for ci in range(len(combinaisons)):
					translations.add(translTable[translTable_id][combinaisons[ci]])
				protein += '/'.join(translations)

	return protein


def polymorphicCodonCombinaisons(codon) :
	"""Returns all the possible amino acids encoded by codon"""
	return getSequenceCombinaisons(codon, 0)

def getSequenceCombinaisons(polymorphipolymorphicDnaSeqSeq, pos = 0) :
	"""Takes a dna sequence with polymorphismes and returns all the possible sequences that it can yield"""

	if type(polymorphipolymorphicDnaSeqSeq) is not types.ListType :
		seq = list(polymorphipolymorphicDnaSeqSeq)
	else :
		seq = polymorphipolymorphicDnaSeqSeq

	if pos >= len(seq) :
		return [''.join(seq)]

	variants = []
	if seq[pos] in polymorphicNucleotides :
		chars = decodePolymorphicNucleotide(seq[pos])
	else :
		chars = seq[pos]#.split('/')

	for c in chars :
		rseq = copy.copy(seq)
		rseq[pos] = c
		variants.extend(getSequenceCombinaisons(rseq, pos + 1))

	return variants

def decodePolymorphicNucleotide(nuc) :
	"""the opposite of encodePolymorphicNucleotide, from 'R' to ['A', 'G']"""
	if nuc in polymorphicNucleotides :
		return polymorphicNucleotides[nuc]

	if nuc in nucleotides :
		return nuc

	raise ValueError('nuc: %s, is not a valid nucleotide' % nuc)