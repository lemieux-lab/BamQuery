from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField, BooleanField, SelectField
from wtforms.validators import DataRequired, Length, Optional
from flask_wtf.file import FileField, FileAllowed, FileRequired


class MyFloatField(FloatField):
	def process_formdata(self, valuelist):
		if valuelist:
			try:
				self.data = float(valuelist[0].replace(',', '.'))
			except ValueError:
				self.data = None
				raise ValueError(self.gettext('Not a valid float value'))

# https://getbootstrap.com/docs/5.0/forms/form-control/
class BamQuery_search(FlaskForm):

	name_query = StringField('Name Query :', 
		validators=[DataRequired(), Length(min=2, max=20)])

	peptides = FileField('Upload peptides.tsv file', 
		validators=[FileRequired(), FileAllowed(['tsv'])])

	strandedness = BooleanField('Strandedness :', default=True)

	th_out = MyFloatField('Threshold :', 
		validators=[Optional()],default=8.55)

	genome_version = SelectField('Genome Version :', choices=['v26_88', 'v33_99', 'v38_104']) 
	data_base_snp = SelectField('dbSNP Release :', choices=['dbSNP_149', 'dbSNP_151', 'dbSNP_155' , 'No_db']) 

	submit = SubmitField('Submit')


class Retrieve_results(FlaskForm):

	name_query_to_retrieve = StringField('Name Query :', 
		validators=[DataRequired(), Length(min=2, max=20)])

	submit_name_query = SubmitField('Submit')

