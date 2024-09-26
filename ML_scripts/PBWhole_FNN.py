# Parameters
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Feedforward neural networks for pseudobulk data')

# Add parameters to the parser
parser.add_argument('--inPath', type=str, help='Folder with the input data')
parser.add_argument('--resultsPath', type=str, help='Folder to save the generated reports')
parser.add_argument('--varColumn', type=str, default='Condition', help='Column in Phenotype.tsv that contains the analyzed variable')
parser.add_argument('--sampleColumn', type=str, default='Sample', help='Column in the metadata of clusters with the sample name')
parser.add_argument('--KOuter', type=int, default=5, help='Folds for the outer cross-validation')
parser.add_argument('--KInner', type=int, default=4, help='Folds for the inner cross-validation')
parser.add_argument('--cores', type=int, default=4, help='Number of threads')

# Parse the arguments
args = parser.parse_args()

# Access the parameter values
inPath = args.inPath
resultsPath = args.resultsPath
varColumn = args.varColumn
sampleColumn = args.sampleColumn
KOuter = args.KOuter
KInner = args.KInner
cores = args.cores

# Import libraries
import sys
import pandas as pd
import numpy as np
import os
import sklearn.metrics as metr
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import make_scorer
import glob
import itertools
import statistics
from sklearn.model_selection import cross_validate
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
import torch
from skorch import NeuralNetClassifier
from torch import nn


# Create the output folders if necessary
if not os.path.exists(resultsPath):
    os.mkdir(resultsPath)

test_Results = {}
testPredictions = {}
validationMCCs = {}
trainedModels = {}
resultsPath += "/"

metadataSamples = pd.read_table(inPath + '/Phenodata.tsv')



## Assign categorical labels to numbers
labels = sorted(list(set(metadataSamples[varColumn])))
labelsDict = {}
x = 0

for label in labels:
	labelsDict[label] = x
	x += 1

metadataSamples["LabelInt"] = metadataSamples[varColumn].map(labelsDict)



# Read files
expression = pd.read_table(inPath + 'pseudobulk_whole.tsv').transpose()
metadata = metadataSamples.loc[expression.index.tolist()]

# Outer stratified cross-validation
outer_cv = StratifiedKFold(n_splits=KOuter, shuffle=True, random_state=123)

# Initialize lists to store results
best_params_dict = {}
MCCs = []

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



nGenes = expression.shape[1]
layer1 = round(nGenes/2)
layer2 = round(layer1/2)
layer3 = round(layer2/2)
layer4 = round(layer3/4)
outNeurons = len(labels)
            
net = NeuralNetClassifier(
    NeuralNetwork(layer1, layer2, layer3, layer4, outNeurons = outNeurons, nGenes = nGenes),
    max_epochs=10,
    criterion=nn.CrossEntropyLoss,
    lr=0.1,
    batch_size=10,
    # Shuffle training data on each epoch
    iterator_train__shuffle=True,
    device='cpu',
)
# X = np.array(expression).astype(np.float32)
# y = np.array(metadata["LabelInt"])
# net.fit(X, y)
# y_pred = net.predict(X)

param_grid = {
    'lr': [0.1, 0.01, 0.001, 0.0001, 0.00001], 
    'batch_size': [10, 20, 30, 40, 50, 60],
    'max_epochs': range(10, 310, 20), 
}

# Outer stratified cross-validation loop
for foldOut, (train_index, test_index) in enumerate(outer_cv.split(expression, metadata["LabelInt"].tolist())):
    X_train, X_test = np.array(expression.iloc[train_index].astype(np.float32)), np.array(expression.iloc[test_index].astype(np.float32))
    y_train, y_test = metadata.iloc[train_index,:]["LabelInt"], metadata.iloc[test_index,:]["LabelInt"]

    # Inner cross-validation
    inner_cv = StratifiedKFold(n_splits=KInner, shuffle=True, random_state=123)
    
    # Randomized search for hyperparameter tuning
    random_search = RandomizedSearchCV(
        net, param_grid, n_iter=100, scoring=make_scorer(matthews_corrcoef), cv=inner_cv, n_jobs=cores,
        random_state=123
    )
    
    # Fit the random search on the inner training data
    _ = random_search.fit(X_train, y_train)
    
    # Obtain the best parameters and score
    best_params = random_search.best_params_
    best_params_dict[foldOut] = best_params
    MCCs.append(random_search.best_score_)
    
    # Build a new model using the best parameters and evaluate on the outer test data
    best_model = NeuralNetClassifier(
        NeuralNetwork(layer1, layer2, layer3, layer4, outNeurons = outNeurons, nGenes = nGenes),
        max_epochs=best_params['max_epochs'],
        criterion=nn.CrossEntropyLoss,
        lr=best_params['lr'],
        batch_size=best_params['batch_size'],
        # Shuffle training data on each epoch
        train_split=None,
        device='cpu'
    )
    _ = best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_test)
    for sample in range(len(y_pred)):
        sampleName = y_test.index[sample]
        testPredictions[sampleName] = y_pred[sample]
    

# Train and save the model with the whole data for external validation
# Hyperparameters are the ones from the fold with the best performance
bestFold = MCCs.index(max(MCCs))
best_params = best_params_dict[bestFold]
modelCluster = best_model = NeuralNetClassifier(
    NeuralNetwork(layer1, layer2, layer3, layer4, outNeurons = outNeurons, nGenes = nGenes),
    max_epochs=best_params['max_epochs'],
    criterion=nn.CrossEntropyLoss,
    lr=best_params['lr'],
    batch_size=best_params['batch_size'],
    # Shuffle training data on each epoch
    train_split=None,
    device='cpu',
)
_ = modelCluster.fit(np.array(expression.astype(np.float32)), metadata["LabelInt"])

trainedModels['Model'] = _




# Estimate performance
    
labelsReal = metadataSamples["LabelInt"].loc[list(testPredictions.keys())]
x = list(labelsReal)
y = list(testPredictions.values())
test_Results = {'accuracy': metr.accuracy_score(x, y),
               'precision': metr.precision_score(x, y),
               'recall': metr.recall_score(x, y),
               'f1': metr.f1_score(x, y),
               'MCC': matthews_corrcoef(x, y)}

######################
### Export results ###
######################
results_pd = pd.DataFrame.from_dict(test_Results, orient = "index")
results_pd.to_csv(resultsPath + "testResults.tsv", sep="\t")

torch.save(trainedModels, resultsPath + "models.pt")
