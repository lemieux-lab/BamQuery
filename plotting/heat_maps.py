import os, logging, time, subprocess, pickle, multiprocessing, os, _thread, csv, collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style

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
	plt.savefig(path_to_output_folder+'plots/'+name_exp+'_'+name+'.pdf', format='pdf', bbox_inches='tight', pad_inches=0, orientation='landscape')
	plt.show()