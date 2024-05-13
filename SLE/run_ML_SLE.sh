#######################
### Pseudobulk data ###
#######################

Rscript ML_scripts/Pseudobulk.R --inputPath SLE/data/SLE_Science/ --sampleColumn ind_cov
Rscript ML_scripts/Pseudobulk.R --inputPath SLE/data/SLE_pediatrics/ --sampleColumn ind_cov

###################
### Run models  ###
###################

algorithms=("DT" "KNN" "LDA" "LR" "NB" "RF" "SVM")

for algorithm in "${algorithms[@]}"; do
    # Training and internal validation with pseudobulk
    python 'ML_scripts/PBWhole_'$algorithm'.py' --inPath SLE/data/SLE_Science/ --sampleColumn ind_cov  --resultsPath 'SLE/results_SLE/'$algorithm'_Whole/' --varColumn Status  --KOuter 5 --KInner 4 --cores 128
    
    # External validation
    python ML_scripts/PBWhole_pretrained.py --inPath SLE/data/SLE_pediatrics/ --modelFile 'SLE/results_SLE/'$algorithm'_Whole/models.pt' --sampleColumn ind_cov  --outFile 'SLE/results_SLE/'$algorithm'_Whole/pediatrics_prediction.tsv'
done
