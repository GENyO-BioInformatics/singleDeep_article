# Parameters
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Random Forest models for pseudobulk data')

# Add parameters to the parser
parser.add_argument('--inPath', type=str, help='Folder with the input data')
parser.add_argument('--modelFile', type=str, help='Saved model')
parser.add_argument('--outFile', type=str, default='prediction_results.tsv', help='Name of the output table')
parser.add_argument('--sampleColumn', type=str, default='Sample', help='Column in the metadata of clusters with the sample name')

# Parse the arguments
args = parser.parse_args()

# Access the parameter values
inPath = args.inPath
modelFile = args.modelFile
outFile = args.outFile
sampleColumn = args.sampleColumn

# Import libraries
import pandas as pd
import glob
import itertools
import torch

labelsPredicted = {}

metadataSamples = pd.read_table(inPath + '/Phenodata.tsv')

training_dict = torch.load(modelFile)


expression = pd.read_table(inPath + '/pseudobulk_whole.tsv').transpose()
prediction = training_dict['Model'].predict(expression)
for sample in range(len(prediction)):
    sampleName = expression.index[sample]
    labelsPredicted[sampleName] = prediction[sample]

######################
### Export results ###
######################
results_pd = pd.DataFrame.from_dict(labelsPredicted, orient = "index")
results_pd.columns = ['label_predicted']
results_pd.to_csv(outFile, sep="\t")
