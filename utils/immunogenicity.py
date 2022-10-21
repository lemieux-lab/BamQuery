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

		rphm_norm = pd.read_csv(data)

		data_3 = []

		for i, row in rphm_norm.iterrows():
			peptide_type = row['Peptide Type']
			peptide = row['Peptide']
			samples_mtec = ['star_062015', 'star_102015', 'star_S5', 'star_S9', 'star_S10', 'star_S11', 'star_S16_mTECs', 'star_S15_mTECs', 'star_ERR2882510', 'star_ERR2882513', 'star_ERR2882514']
			values_mtec = row[samples_mtec].values
			mean_mtec = np.mean(values_mtec)
			
			samples_blood = ['star_SRR7300571', 'star_SRR7300573', 'star_SRR7300593', 'star_SRR7300595', 'star_SRR7300613', 'star_SRR7300651', 'star_SRR7300653', 'star_SRR3084021', 'star_SRR3084022', 'star_SRR3084023', 'star_SRR3084024', 'star_SRR3084025', 'star_SRR3084026', 'star_SRR3084027', 'star_SRR3084028', 'star_SRR3084029', 'star_SRR3084030', 'star_SRR3084031', 'star_SRR3084032']
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





	