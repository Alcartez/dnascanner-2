from flask_wtf import FlaskForm
from wtforms import TextAreaField, EmailField, BooleanField, SubmitField, IntegerField
from flask_wtf.file import FileField
import pandas as pd

class FastaForm(FlaskForm):
    fasta_text = TextAreaField('Enter Nucleotide Sequence in Fasta Field')
    fasta_file = FileField('Choose File')
    email = EmailField('Email')
    windowWidth = IntegerField("Enter Window Width")
    submit = SubmitField("Submit")

param_file_2 = "Parameter_Files/Parameter_Files_Trinucleotide - Sheet1.csv"
param_file_1 = "Parameter_Files/Parameter_Sheet_Dinucleotide - Sheet1.csv"

parameters = {
    'Dinucleotide' : pd.read_csv(param_file_1).columns.values.tolist()[1:],
    'Trinucleotide' : pd.read_csv(param_file_2).columns.values.tolist()[1:]
}