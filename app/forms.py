from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField, BooleanField, SelectField, FieldList
from wtforms.validators import DataRequired, Length, Optional
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import TextAreaField
from wtforms.widgets import TextArea
from wtforms.validators import ValidationError
import re

aa = ['A', 'G', 'P', 'T', 'V', 'C', 'D', 'E', 'F', 'H', 'K', 'N', 'Q', 'Y', 'I', 'L', 'R', 'S', 'M', 'W']

class MyFloatField(FloatField):
	def process_formdata(self, valuelist):
		if valuelist:
			try:
				self.data = float(valuelist[0].replace(',', '.'))
			except ValueError:
				self.data = None
				raise ValueError(self.gettext('Not a valid float value'))

#https://wtforms.readthedocs.io/en/stable/validators/
def validate_peptides_list(form, field):
		peptides = field.data.split('\n')
		m = re.compile(r'[AGPTVCDEFHKNQYILRSMW]')
		for pep in peptides:
			if pep.strip().isdigit():
				raise ValidationError('Number found in peptide sequence.')
			if len(pep.strip()) > 11:
				print (pep.strip(), len(pep.strip()))
				raise ValidationError('BamQuery analyze MAPs ranging in length from 8 to 11 amino acids (aa).')
			if not m.match(pep.strip()):
				raise ValidationError('Sequence of peptide does not contain valid amino acids (aa).')
		if len(peptides)>100:
			raise ValidationError('Number of peptides allowed per run = 100.')
		

# https://getbootstrap.com/docs/5.0/forms/form-control/
class BamQuery_search(FlaskForm):

	name_query = StringField('Name Query :', 
		validators=[DataRequired(message="This field is required."), Length(min=2, max=20)])

	#peptides = FileField('Upload peptides.tsv file', validators=[FileRequired(), FileAllowed(['tsv'])])

	#peptides_list = StringField('Peptides List', widget=TextArea(), validators=[DataRequired(), Length(min=8, max=11)])

	peptides_list = TextAreaField('Peptides List', validators=[DataRequired(), validate_peptides_list])


	th_out = MyFloatField('Threshold :', 
		validators=[Optional()],default=8.55)

	genome_version = SelectField('Genome Version :', choices=['v26_88', 'v33_99', 'v38_104']) 
	data_base_snp = SelectField('dbSNP Release :', choices=['dbSNP_149', 'dbSNP_151', 'dbSNP_155' , 'None']) 

	submit = SubmitField('Submit')


class Retrieve_results(FlaskForm):

	name_query_to_retrieve = StringField('Enter the name of the query :', 
		validators=[DataRequired(), Length(min=2, max=20)])

	submit_name_query = SubmitField('Retrieve')

