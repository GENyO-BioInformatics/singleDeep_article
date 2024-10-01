# Parameters
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='KNN models for pseudobulk data')

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
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_validate
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
import torch

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

param_grid = {
    'n_neighbors': [3, 5, 7, 9],
    'weights': ['uniform', 'distance'],
    'algorithm': ['auto', 'ball_tree', 'kd_tree', 'brute'],
    'p': [1, 2],
}

# Read files
expression = pd.read_table(inPath + 'pseudobulk_whole.tsv').transpose()
metadata = metadataSamples.loc[expression.index.tolist()]

# Outer stratified cross-validation
outer_cv = StratifiedKFold(n_splits=KOuter, shuffle=True, random_state=123)

# Initialize lists to store results
best_params_dict = {}
MCCs = []

# Outer stratified cross-validation loop
for foldOut, (train_index, test_index) in enumerate(outer_cv.split(expression, metadata["LabelInt"].tolist())):
    X_train, X_test = expression.iloc[train_index], expression.iloc[test_index]
    y_train, y_test = metadata.iloc[train_index,:]["LabelInt"], metadata.iloc[test_index,:]["LabelInt"]

    # Inner cross-validation
    inner_cv = StratifiedKFold(n_splits=KInner, shuffle=True, random_state=123)
    
    # Randomized search for hyperparameter tuning
    knn_model = KNeighborsClassifier()
    random_search = RandomizedSearchCV(
        knn_model, param_grid, n_iter=100, scoring=make_scorer(matthews_corrcoef), cv=inner_cv, n_jobs=cores,
        random_state=123
    )
    
    # Fit the random search on the inner training data
    _ = random_search.fit(X_train, y_train)
    
    # Obtain the best parameters and score
    best_params = random_search.best_params_
    best_params_dict[foldOut] = best_params
    MCCs.append(random_search.best_score_)
    
    # Build a new model using the best parameters and evaluate on the outer test data
    best_model = KNeighborsClassifier(**best_params)
    _ = best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_test)
    for sample in range(len(y_pred)):
        sampleName = y_test.index[sample]
        testPredictions[sampleName] = y_pred[sample]
    

# Train and save the model with the whole data for external validation
# Hyperparameters are the ones from the fold with the best performance
bestFold = MCCs.index(max(MCCs))
best_params = best_params_dict[bestFold]
modelCluster = KNeighborsClassifier(**best_params)
_ = modelCluster.fit(expression, metadata["LabelInt"].tolist())

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
