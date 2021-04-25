import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style
import numpy as np
import random
import scipy.stats
from matplotlib import cm


violetes = sns.cubehelix_palette(7)
blues_palette = sns.color_palette("Blues")
navy_colors = sns.light_palette("navy", reverse=False)[1:]
green_palette = sns.light_palette("green")[1:]

assignation_colors = {'Protein-coding genes': blues_palette[0],
						'Protein-coding transcripts': blues_palette[0],
						'5UTR':  blues_palette[1],
						'3UTR':  blues_palette[2],
						'In_frame':  blues_palette[3],
						'Frameshift':  blues_palette[4],
						'protein_coding':  blues_palette[0],
						'Junctions': blues_palette[5],
						'Other coding regions': blues_palette[5],

						'Non-coding genes': violetes[0],
						'Non-coding transcripts': violetes[0],
						'processed_transcript': violetes[1],
						'nonsense_mediated_decay': violetes[2],
						'antisense': violetes[3],
						'Exons': violetes[4],
						'lincRNA': violetes[5],
						'Other non-coding regions':violetes[6],

						'Intergenic Region': green_palette[0],
						'Intergenic': green_palette[1],
						'Intronic Region': navy_colors[1],
						'Introns':navy_colors[2]}



def get_heat_map(df, path_to_output_folder, name_exp, name, norm, th_out = 8.55):
	peptides_total = len(df.index)
	bam_files = len(df.columns)

	width = 10
	heigth = 10
	fontsize = 12
	annot = True

	if peptides_total > 60 :
		heigth = 20
		fontsize = 8
		annot = False

	if bam_files > 30 :
		width = 20
		fontsize = 8
		annot = False

	fig = plt.figure(figsize=(width, heigth), dpi=300)
	ax = fig.add_subplot(111)
	ax.grid(False)
	
	if 'norm' not in name :
		ax = sns.heatmap(df, cmap="Blues",linewidths=1, linecolor='white',xticklabels = 1, yticklabels = 1, annot=annot, fmt='g') #cmap="YlGnBu"
	else:
		ax = sns.heatmap(df, cmap="Blues",linewidths=1, linecolor='white',xticklabels = 1, yticklabels = 1, annot=annot, fmt='.1f') #cmap="YlGnBu"
	plt.yticks(rotation=0)
	plt.tight_layout()
	plt.savefig(path_to_output_folder+'plots/heat_maps/'+name_exp+name+'.pdf', format='pdf', bbox_inches='tight', pad_inches=0, orientation='landscape')
	plt.show()

	if norm:
		script_R_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/average_tissues_mode.R'
		command = 'Rscript '+script_R_path+' '+path_to_output_folder+'res/AUX_files/processed/'+' '+path_to_output_folder+'plots/heat_maps/ '+str(th_out)+' '+name_exp+name
		subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, close_fds=True)


def plot_pie(title, outer_labels, intra_labels, intra_sizes, outer_sizes, path_to_output_folder, name_exp, name, fontsize=12):
	
	fig, ax = plt.subplots(figsize=(11, 7), dpi=300) 

	title=title+'\n'

	colors_outer_labels = []
	colors_intra_labels = []

	for label in outer_labels:
		colors_outer_labels.append(assignation_colors[label])

	for label in intra_labels:
		colors_intra_labels.append(assignation_colors[label])

	# Plot
	plt.title(title, fontsize = fontsize+2)

	wedges1, texts1 = plt.pie(outer_sizes, colors=colors_outer_labels, 
							startangle=180, pctdistance=0.9, shadow=False, 
							textprops={'fontsize': fontsize, 'fontweight' : 'normal'},
							wedgeprops={"edgecolor":"w",'linewidth': 5, 'antialiased': True})

	wedges, texts, autotexts = plt.pie(intra_sizes, colors=colors_intra_labels, 
							radius=0.8,startangle=180, 
							textprops={'fontsize': fontsize, 'fontweight' : 'normal'}, 
							shadow=True, autopct='%1.0f%%', pctdistance=0.8,
							wedgeprops={"edgecolor":"w",'linewidth': 1, 'antialiased': True})
	
	#plt.legend(wedges1, outer_labels, loc="best", bbox_to_anchor=(0.7,-0.05, 0.5, 0.5), fontsize = fontsize, frameon=True)

	plt.legend(wedges, intra_labels, loc="best", bbox_to_anchor=(0.7,-0.05, 0.5, 0.8), fontsize = fontsize, frameon=True)
	
	bbox_props = dict(boxstyle="square,pad=0.2", color='black', fc="w", ec="k", lw=0.2)
	kw = dict(arrowprops=dict(arrowstyle="-", color='black', linewidth=2),zorder=2, va="center")
	plt.axis('equal')

	#plt.savefig(path_to_output_folder+'plots/biotypes/'+name_exp+'_'+name+'.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.5)
	#plt.show()

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
			plt.savefig(path_to_output_folder+name_exp+'_'+name+'.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.5)
			plt.show()
	except :
		print ('It was not possible to save the figure for ', title )


def plot_pie_ere(title, outer_labels, outer_sizes, path_to_output_folder, name_exp, name, fontsize=12):
	
	number = len(outer_labels)
	#a = np.random.random(number)
	#cs = cm.Set1(np.arange(number)/(number * 1.0))
	
	import matplotlib.colors as mcolors
	colors = random.choices(list(mcolors.CSS4_COLORS.values()),k = number)

	fig, ax = plt.subplots(figsize=(11, 7), dpi=300) 

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
	plt.show()


def correlation(path_to_output_folder, name_exp, dataframe):

	slope, intercept, r, p, stderr = scipy.stats.linregress(dataframe['Read count RNA'], dataframe['Read count Ribo'])
	rho, p = scipy.stats.pearsonr(dataframe['Read count RNA'], dataframe['Read count Ribo'])
	x = dataframe['Read count RNA']

	fig, ax = plt.subplots(figsize=(10,6), dpi=300)
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
	plt.show()



	



	
