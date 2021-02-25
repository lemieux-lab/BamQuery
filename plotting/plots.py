import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style
import numpy as np
import random

violetes = sns.cubehelix_palette(7)
blues_palette = sns.color_palette("Blues")
navy_colors = sns.light_palette("navy", reverse=False)[1:]
green_palette = sns.light_palette("green")[1:]
#others = 'retained_intron', 'IG_V_gene', 'TEC', 'bidirectional_promoter_lncRNA', 'transcribed_unitary_pseudogene', 
#'transcribed_unprocessed_pseudogene', 'sense_overlapping','processed_pseudogene', 'unprocessed_pseudogene'

assignation_colors = {'Protein-coding genes': blues_palette[0],
						'5UTR':  blues_palette[1],
						'3UTR':  blues_palette[2],
						'In_frame':  blues_palette[3],
						'Frameshift':  blues_palette[4],
						'protein_coding':  blues_palette[0],
						'CDS': blues_palette[5],

						'Non-coding genes': violetes[0],
						'processed_transcript': violetes[1],
						'nonsense_mediated_decay': violetes[2],
						'antisense': violetes[3],
						'Exons': violetes[4],
						'lincRNA': violetes[5],
						'Others':violetes[6],

						'Intergenic': green_palette[0],
						'Intronic': navy_colors[1],
						'Introns':navy_colors[1]}



def get_heat_map(df, path_to_output_folder, name_exp, name):
	fontsize = 12
	fig = plt.figure(figsize=(10, 6), dpi=300)
	ax = fig.add_subplot(111)
	ax.grid(False)
	if 'norm' not in name :
		ax = sns.heatmap(df, cmap="Blues",linewidths=1, linecolor='white',xticklabels = 1, yticklabels = 1, annot=True, fmt='g') #cmap="YlGnBu"
	else:
		ax = sns.heatmap(df, cmap="Blues",linewidths=1, linecolor='white',xticklabels = 1, yticklabels = 1, annot=True, fmt='.1f') #cmap="YlGnBu"
	plt.yticks(rotation=0)
	plt.tight_layout()
	plt.savefig(path_to_output_folder+'plots/heat_maps/'+name_exp+'_'+name+'.pdf', format='pdf', bbox_inches='tight', pad_inches=0, orientation='landscape')
	plt.show()


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
	
	plt.legend(wedges1, outer_labels, loc="best", bbox_to_anchor=(0.7,-0.05, 0.5, 0.5), fontsize = fontsize, frameon=True)
	
	bbox_props = dict(boxstyle="square,pad=0.2", color='black', fc="w", ec="k", lw=0.2)
	kw = dict(arrowprops=dict(arrowstyle="-", color='black', linewidth=2),zorder=2, va="center")
	plt.axis('equal')

	try:
		for i, p in enumerate(wedges):
			ang = (p.theta2 - p.theta1)/2. + p.theta1 
			if ang == 360:
				ang = 359

			y = np.sin(np.deg2rad(ang))
			x = np.cos(np.deg2rad(ang))

			horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
			connectionstyle = "angle,angleA=0,angleB={}".format(ang)
			kw["arrowprops"].update({"connectionstyle": connectionstyle})
			
			ax.annotate(intra_labels[i], xy=(x,y), xytext=(1.2*np.sign(x), 1.2*y), horizontalalignment=horizontalalignment, fontsize = fontsize, fontweight = 'normal', **kw)
			plt.savefig(path_to_output_folder+'plots/biotypes/'+name_exp+'_'+name+'.pdf', orientation='landscape', format='pdf', bbox_inches='tight', pad_inches=0.5)
			plt.show()
	except :
		print ('It was not possible to save the figure for ', title )

	
