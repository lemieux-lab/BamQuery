import warnings
warnings.filterwarnings("ignore")
import math
import numpy as np
import os
import pandas as pd
import pickle


path_to_lib = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/lib/'


class Immunogenicity:

	def __init__(self, path_to_output_folder, name_exp):
		self.path_to_output_folder = path_to_output_folder
		self.name_exp = name_exp


	def get_predictions(self):

		# Make prediction 
		data = self.path_to_output_folder+'/res/temps_files/'+self.name_exp+'_rna_norm.csv'
		model_path = path_to_lib+'logistic_regression_best_model_cv_aug_31.sav'
		regression_model = pickle.load(open(model_path, 'rb'))

		info_bam_files_tissues = self.path_to_output_folder+'/res/info_bam_files_tissues.csv'
		info_bam_files_tissues = pd.read_csv(info_bam_files_tissues) 


		rphm_norm = pd.read_csv(data)

		data_3 = []

		for i, row in rphm_norm.iterrows():
			peptide_type = row['Peptide Type']
			peptide = row['Peptide']
			samples_mtec = info_bam_files_tissues.loc[(info_bam_files_tissues['Sample category'] == 'mTEC')]['sample_ids'].values[0].split(' ')
			values_mtec = row[samples_mtec].values
			mean_mtec = np.mean(values_mtec)
			
			samples_blood = info_bam_files_tissues.loc[(info_bam_files_tissues['Sample category'] == 'Blood_DC')]['sample_ids'].values[0].split(' ')
			values_blood = row[samples_blood].values
			mean_blood = np.mean(values_blood)
			
			data_3.append([peptide_type, peptide, mean_mtec, mean_blood])

		df_summary_3 = pd.DataFrame(data_3, columns=['Peptide Type','Peptide','mTEC_RPHM', 'DC_RPHM'])
		X = df_summary_3[['mTEC_RPHM', 'DC_RPHM']]
		ypred = regression_model.predict(X)
		y_pred_proba = regression_model.predict_proba(X)[::,1]

		df_summary_3['Probability'] = y_pred_proba
		df_summary_3['Immunogenecity'] = ypred
		df_summary_3['Immunogenecity'] = df_summary_3['Immunogenecity'].map({1.0: 'Immunogenic', 0.0: 'Non_immunogenic'})
		
		path_output_prediction = self.path_to_output_folder +'/res/immunogenicity_prediction.csv'
		df_summary_3.to_csv(path_output_prediction, header=True, index=False)





	