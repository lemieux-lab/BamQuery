import plotting.plots as plots
import gc

def draw_biotypes(biotypes_peptides, path_to_save, global_, samples, name_exp):

	others_non_coding = ['retained_intron', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene']
	others_protein_coding = ['IG_V_gene', 'TEC', 'Exons']
	others_ere = ['DNA','RC', 'RNA','Satellite','Simple_repeat','Unknown', 'Low_complexity', 'rRNA','scRNA','snRNA','srpRNA','tRNA']
	

	organisation_labels = {'Protein-coding Regions':['5UTR', '3UTR', 'In_frame', 'Frameshift', 'protein_coding', 'CDS', 'Junctions', 'Other coding regions'], 
						'Non-coding RNAs':['processed_transcript', 'nonsense_mediated_decay', 'antisense', 'Non_coding Exons', 'Non_coding Junctions', 'lincRNA', 'Other non_coding regions'], 
						'Intergenic Regions':['Intergenic'], 
						'Intronic Regions':['Introns'],
						'EREs':['LINE', 'LTR','Retroposon','SINE', 'Other EREs']}

	def plot_biotype(biotypes, name):
		title = name
		labels_in_type_peptide = {'Protein-coding genes':{}, 'Non-coding genes': {}, 'Protein-coding Regions':{}, 'Non-coding RNAs': {}, 'Protein-coding transcripts':{}, 'Non-coding transcripts': {}, 'Intergenic Regions':{}, 'Intronic Regions':{}, 'EREs':{}}
		outer_labels = []
		outer_sizes = []
		intra_labels = []
		intra_sizes = []

		for biotype, total_biotype in biotypes.items():

			in_ = False
			for type_, types in organisation_labels.items():
				if biotype in types:
					labels_in_type_peptide[type_][biotype] = total_biotype 
					in_ = True
					break

			if not in_:
				
				if biotype in others_non_coding:
					type_ = 'Non-coding Regions'
					try:
						labels_in_type_peptide[type_]['Other non-coding regions'] += total_biotype 
					except KeyError:
						labels_in_type_peptide[type_]['Other non-coding regions'] = total_biotype

				elif biotype in others_protein_coding:
					type_ = 'Protein-coding Regions'
					try:
						labels_in_type_peptide[type_]['Other coding regions'] += total_biotype
					except KeyError:
						labels_in_type_peptide[type_]['Other coding regions'] = total_biotype

				elif biotype in others_ere:
					type_ = 'EREs'
					try:
						labels_in_type_peptide[type_]['Other EREs'] += total_biotype
					except KeyError:
						labels_in_type_peptide[type_]['Other EREs'] = total_biotype

				else:
					if '-' in biotype:
						if 'Non_coding' in biotype :
							try:
								labels_in_type_peptide['Non-coding RNAs']['Non_coding Junctions'] += total_biotype 
							except KeyError:
								labels_in_type_peptide['Non-coding RNAs']['Non_coding Junctions'] = total_biotype
						else:
							try:
								labels_in_type_peptide['Protein-coding Regions']['Junctions'] += total_biotype 
							except KeyError:
								labels_in_type_peptide['Protein-coding Regions']['Junctions'] = total_biotype
					else:
						print ('Problem with assignation biotype : ', biotype, name)
						
					
		for outer_label_type, intra_labels_type in organisation_labels.items():
			values_in = 0
			for intra_label in intra_labels_type:
				try:
					value = labels_in_type_peptide[outer_label_type][intra_label]
					intra_labels.append(intra_label)
					intra_sizes.append(value)
					values_in += value
				except:
					pass
			if values_in > 0:
				outer_labels.append(outer_label_type)
				outer_sizes.append(values_in)
			
		plots.plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_save, name_exp, name)
		gc.collect()

	if not samples:
		if not global_ :
			
			for type_peptide, biotypes in biotypes_peptides.items():
				plot_biotype(biotypes, type_peptide)
		else:
			plot_biotype(biotypes_peptides, 'All_peptides')
	else:
		if not global_  :
			for key_group, key_peptides_dic in biotypes_peptides.items():
				for type_peptide, biotypes in key_peptides_dic.items():
					plot_biotype(biotypes, type_peptide+'_'+key_group)
		else:
			for key_group, biotypes in biotypes_peptides.items():
				plot_biotype(biotypes, 'All_peptides_'+key_group)

def draw_correlation(data, name_exp, path_to_output_folder):

	#[type_peptide, peptide, alignment, strand, 'No Annotation', gene_level_biotype, transcript_level_biotype, genomic_position_biotype, count_rna, count_ribo]
	counts_ribo = []
	counts_rna = []
	peptide_type = []
	peptides = []

	for peptide_type_, peptides_ in data.items():

		for peptide, counts in peptides_.items():
			counts_ribo.append(counts[1])
			counts_rna.append(counts[0])
			peptides.append(peptide)
		peptide_type.extend([peptide_type_]*len(peptides_))
		
	data = {'Type Peptide': peptide_type, 
			'Peptides': peptides,
    		'Read count Ribo': counts_ribo,
    		'Read count RNA': counts_rna} 
    		
	df_correlation = pd.DataFrame(data)
	plots.correlation(path_to_output_folder, name_exp, df_correlation)

