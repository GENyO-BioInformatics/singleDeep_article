# Parameters
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Pretrained models for pseudobulk data')

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
import numpy as np
import glob
import itertools
import torch
from torch import nn

# Network class
class NeuralNetwork(nn.Module):
    def __init__(self, Hs1, Hs2, Hs3, Hs4, outNeurons, nGenes):
        super(NeuralNetwork, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(nGenes, Hs1),
            nn.ReLU(), # ReLU activation function
            nn.Linear(Hs1, Hs2),
            nn.ReLU(),
            nn.Linear(Hs2, Hs3),
            nn.ReLU(),
            nn.Linear(Hs3, Hs4),
            nn.ReLU(),
            nn.Linear(Hs4, outNeurons),
        )

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits


labelsPredicted = {}

metadataSamples = pd.read_table(inPath + '/Phenodata.tsv')

training_dict = torch.load(modelFile)

expression = pd.read_table(inPath + '/pseudobulk_whole.tsv').transpose()

if 'FNN' in modelFile:
    expression = np.array(expression.astype(np.float32))
    prediction = training_dict['Model'].predict(expression)
    expression = pd.read_table(inPath + '/pseudobulk_whole.tsv').transpose()
else:
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
