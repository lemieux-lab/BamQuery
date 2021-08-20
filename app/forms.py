from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField, BooleanField
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

class BamQuery_search(FlaskForm):

	name_query = StringField('Name Query :', 
		validators=[DataRequired(), Length(min=2, max=20)])

	peptides = FileField('Upload peptides.tsv file', 
		validators=[FileRequired(), FileAllowed(['tsv'])])

	strandedness = BooleanField('Strandedness')

	th_out = MyFloatField('Threshold :', 
		validators=[Optional()])

	submit = SubmitField('Submit')

