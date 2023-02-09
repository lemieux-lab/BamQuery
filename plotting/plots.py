import os, subprocess, pickle, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats
import matplotlib.colors as mcolors
	

__author__ = "Maria Virginia Ruiz Cuevas"
__email__ = "maria.virginia.ruiz.cuevas@umontreal.ca"

violetes = sns.cubehelix_palette(7)
blues_palette = sns.color_palette("Blues",n_colors = 8)
navy_colors = sns.light_palette("navy", reverse=False)[1:]
green_palette = sns.light_palette("green")[1:]
flare_palette = sns.color_palette("flare", n_colors = 6)
tab10_palette = sns.color_palette("tab10")

assignation_colors = {  'Protein-coding Regions': blues_palette[0],
						'5UTR':  blues_palette[1],
						'3UTR':  blues_palette[2],
						'In_frame':  blues_palette[3],
						'Frameshift':  blues_palette[4],
						'CDS':  blues_palette[5],
						'Junctions': blues_palette[6],
						'Other coding regions': blues_palette[7],

						'Non-coding RNAs': violetes[0],
						'Non_coding Exons': violetes[1],
						'Non_coding Junctions': violetes[2],
						
						'Intergenic Regions': green_palette[0],
						'Intergenic': green_palette[1],

						'Intronic Regions': navy_colors[1],
						'Introns':navy_colors[2],

						'EREs' : flare_palette[0],
						'LINE': flare_palette[1],
						'LTR': flare_palette[2], 
						'SINE': flare_palette[3],
						'Antisense_EREs': flare_palette[4],
						'Other EREs': flare_palette[5] }

def get_heat_map(df, path_to_output_folder, mode, path_temps_file, name_exp, name, norm, th_out = 8.55):

	peptides_total = len(df.index)
	bam_files = len(df.columns)
	path_to_output_folder_aux = path_to_output_folder

	if peptides_total < 400 and bam_files < 200:
		
		data = []
		df.reset_index()
		columns_names_bam_files = list(df.columns)[2:]
		
		with open(path_temps_file+'/bam_files_info_query.dic', 'rb') as handle:
			bam_files_info_query = pickle.load(handle)

		# This is to put for each peptide type all the peptides corresponding to them. It doesn't take into account that
		# the peptide can be assoacited with other peptide type

		# for index, row in df.iterrows():
		# 	peptide_type = row['Peptide Type'].split(';')
		# 	for type_ in peptide_type:
		# 		peptide = row['Peptide']
		# 		aux = [type_, peptide]
		# 		for i, sample in enumerate(columns_names_bam_files):
		# 			aux.append(row[sample])
		# 			aux.append(sample)
		# 			tissue = bam_files_info_query[sample][6]
		# 			aux.append(tissue)
		# 			data.append(aux)
		# 			aux = [type_, peptide]

		for index, row in df.iterrows():
			peptide_type = row['Peptide Type']
			peptide = row['Peptide']
			aux = [peptide_type, peptide]
			for i, sample in enumerate(columns_names_bam_files):
				aux.append(row[sample])
				aux.append(sample)
				tissue = bam_files_info_query[sample][6]
				aux.append(tissue)
				data.append(aux)
				aux = [peptide_type, peptide]
		
		if mode == 'translation':
			path_to_output_folder = path_to_output_folder_aux+'/total_translation_expression_heatmap/'
		else:
			path_to_output_folder = path_to_output_folder_aux+'/total_transcription_expression_heatmap/'

		df_counts = pd.DataFrame(data, columns=['Peptide_Type', 'Peptide', 'value', 'Sample', 'Tissue'])	
		path = path_to_output_folder+name_exp+name+'.csv'
		script_R_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/total_expression.R'
		
		if norm:
			label = 'Log10\(RPHM\+1\)'
		else:
			label = 'Log10\(Reads_count+1\)'
			df_counts['value'] = df_counts['value'].apply( lambda x: np.log10(x+1))
			
		df_counts.to_csv(path, index=False)

		command = 'Rscript '+script_R_path+' '+path+' '+path_to_output_folder+' '+name_exp+name+' '+label+' '+str(th_out)
		subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)

	if norm and peptides_total < 400:

		if mode == 'translation':
			path_to_output_folder = path_to_output_folder_aux+'/average_translation_expression_heatmap/'
		else:
			path_to_output_folder = path_to_output_folder_aux+'/average_transcription_expression_heatmap/'

		script_R_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/average_tissues_mode.R'
		command = 'Rscript '+script_R_path+' '+path_to_output_folder+'/norm_info.csv ' +path_to_output_folder+' '+str(th_out)+' '+name_exp+name+' '+mode
		subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True).wait()


def get_heat_map_coverage(df, path_to_output_folder, path_temps_file, name_exp, name):

	peptides_total = len(df.index)
	bam_files = len(df.columns)
	
	if peptides_total < 400 and bam_files < 200:
		
		data = []
		df.reset_index()
		columns_names_bam_files = list(df.columns)[2:]
		
		with open(path_temps_file+'/bam_files_info_query.dic', 'rb') as handle:
			bam_files_info_query = pickle.load(handle)

		for index, row in df.iterrows():
			peptide_type = row['Peptide Type']
			peptide = row['Peptide']
			aux = [peptide_type, peptide]
			for i, sample in enumerate(columns_names_bam_files):
				aux.append(row[sample])
				aux.append(sample)
				tissue = bam_files_info_query[sample][6]
				aux.append(tissue)
				data.append(aux)
				aux = [peptide_type, peptide]

		df_tpm = pd.DataFrame(data, columns=['Peptide_Type', 'Peptide', 'value', 'Sample', 'Tissue'])	

		path = path_to_output_folder+name_exp+name+'.csv'
		df_tpm.to_csv(path, index=False)

		script_R_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/total_expression.R'
		command = 'Rscript '+script_R_path+' '+path+' '+path_to_output_folder+' '+name_exp+name+" Log10\(TPM\+1\)"
		subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True).wait()


def plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_output_folder, name_exp, name, fontsize=12):
	
	fig, ax = plt.subplots(figsize=(11, 7)) # dpi=300) 

	title=title+'\n'

	colors_outer_labels = []
	colors_intra_labels = []

	for label in outer_labels:
		colors_outer_labels.append(assignation_colors[label])

	for label in intra_labels:
		colors_intra_labels.append(assignation_colors[label])

	# Plot
	plt.title(title, fontsize = fontsize+2)

	wedges1, texts1, autotexts1 = plt.pie(outer_sizes, colors=colors_outer_labels, 
							startangle=180, pctdistance=0.9, shadow=False, 
							textprops={'fontsize': fontsize, 'fontweight' : 'normal'},
							autopct='%1.0f%%', wedgeprops={"edgecolor":"w",'linewidth': 1, 'antialiased': True})

	wedges, texts, autotexts = plt.pie(intra_sizes, colors=colors_intra_labels, 
							radius=0.8, startangle=180, 
							textprops={'fontsize': fontsize, 'fontweight' : 'normal'}, 
							shadow=True, autopct='%1.1f%%', pctdistance=0.8,
							wedgeprops={"edgecolor":"w",'linewidth': 1, 'antialiased': True})
	
	#plt.legend(wedges1, outer_labels, loc="best", bbox_to_anchor=(0.7,-0.05, 0.5, 0.5), fontsize = fontsize, frameon=True)

	plt.legend(wedges, intra_labels, loc="best", bbox_to_anchor=(0.7,-0.05, 0.5, 0.8), fontsize = fontsize, frameon=True)
	
	bbox_props = dict(boxstyle="square,pad=0.2", color='black', fc="w", ec="k", lw=0.2)
	kw = dict(arrowprops=dict(arrowstyle="-", color='black', linewidth=2),zorder=2, va="center")
	plt.axis('equal')

	
	try:
		for i, p in enumerate(wedges1):
			ang = (p.theta2 - p.theta1)/2. + p.theta1 
			if ang == 360:
				ang = 359

			y = np.sin(np.deg2rad(ang))
			x = np.cos(np.deg2rad(ang))

			horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
			connectionstyle = "angle,angleA=0,angleB={}".format(ang)
			kw["arrowprops"].update({"connectionstyle": connectionstyle})
			
			ax.annotate(outer_labels[i], xy=(x,y), xytext=(1.2*np.sign(x), 1.2*y), horizontalalignment=horizontalalignment, fontsize = fontsize, fontweight = 'normal', **kw)
	except :
		print ('It was not possible to save the figure for ', title )

	plt.savefig(path_to_output_folder+name_exp+'_'+name+'.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.5)
	plt.cla()
	plt.close(fig)

def plot_pie_ere(title, outer_labels, outer_sizes, path_to_output_folder, name_exp, name, fontsize=12):
	
	number = len(outer_labels)
	
	colors = random.choices(list(mcolors.CSS4_COLORS.values()),k = number)

	fig, ax = plt.subplots(figsize=(11, 7)) #dpi=300) 

	title=title+'\n'

	plt.title(title, fontsize = fontsize+2)

	wedges1, texts1 = plt.pie(outer_sizes, colors = colors,
							startangle=180, pctdistance=0.9, shadow=False, 
							textprops={'fontsize': fontsize, 'fontweight' : 'normal'},
							wedgeprops={"edgecolor":"black",'linewidth': 1, 'antialiased': True})

	plt.legend(wedges1, outer_labels, loc="best", ncol=2, bbox_to_anchor=(0.7,-0.05, 0, 0.8), fontsize = fontsize, frameon=True)
	
	bbox_props = dict(boxstyle="square,pad=0.2", color='black', fc="w", ec="k", lw=0.2)
	kw = dict(arrowprops=dict(arrowstyle="-", color='black', linewidth=2),zorder=2, va="center")
	plt.axis('equal')

	plt.savefig(path_to_output_folder+name_exp+'_'+name+'.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.5)
	


def correlation(path_to_output_folder, name_exp, dataframe):

	slope, intercept, r, p, stderr = scipy.stats.linregress(dataframe['Read count RNA'], dataframe['Read count Ribo'])
	rho, p = scipy.stats.pearsonr(dataframe['Read count RNA'], dataframe['Read count Ribo'])
	x = dataframe['Read count RNA']

	fig, ax = plt.subplots(figsize=(10,6)) #dpi=300)
	fontsize = 18
	ax.grid(False)
	ax.spines['right'].set_visible(0)
	ax.spines['top'].set_visible(0)
	
	ax = sns.scatterplot(x="Read count RNA", y="Read count Ribo", hue = 'Type Peptide', edgecolor='black', data=dataframe)
	plt.xlabel('Read count RNA', fontsize=fontsize)
	plt.ylabel('Read count Ribo', fontsize=fontsize)

	corr = 'Pearson\'s r = '+str(np.round(r,2))+', Slope = '+str(np.round(slope,2))
	line_1 = ax.plot(x, intercept + slope * x, color='b', label = corr)

	handles, labels = ax.get_legend_handles_labels()
	new_handles = []
	new_labels = []
	
	for i,h in enumerate(handles):
		new_handles.append(h)
		new_labels.append(labels[i])
			
	first_legend = ax.legend(handles=new_handles, labels=new_labels, fancybox=True, fontsize = fontsize-4, frameon=True, shadow=True)
	plt.savefig(path_to_output_folder+'plots/correlation/'+name_exp+'_correlation_ribo_rna.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.3)
	



	



	
