# Explainable deep learning for predicting samples phenotypes from single-cell transcriptomics

Code to reproduce the results from the singleDeep article.

## Prerequisites

Python3 and R should be installed in your system.

Clone the singleDeep repository and install its dependencies:

``` bash
git clone https://github.com/GENyO-BioInformatics/singleDeep.git
pip install -r singleDeep/requirements.txt
```

Create the necessary directories.

```{bash}
mkdir SLE/data/Science figures
```

## Systemic Lupus Erythematosus (SLE) diagnosis

For this use case, a dataset of adult SLE patients and healthy controls (Science dataset) is used to train the models, which are validated with a dataset of pediatric SLE patients and healthy controls (pediatric dataset). Download the Science dataset from <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE174188&format=file&file=GSE174188%5FCLUES1%5Fadjusted%2Eh5ad%2Egz> and store it into the SLE/data/Science/ folder. Download also the pediatrics dataset from <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135779&format=file> and store it into the SLE/data/pediatrics/ folder. The file SLE/data/pediatrics/clinical_info.csv has been extracted from the Supplementary Table 1 of the original article.

It is necessary to read the expression data from the pediatrics dataset and to store it into a h5ad file. For that, use the following command.

```{bash}
Rscript SLE/prepare_pediatrics.R
```

Check that the file SLE/data/pediatrics/pediatrics_raw.h5ad has been created. Then, run the following script to process both datasets, including steps like batch effect correction, regression of mitochondrial genes expression, filtering, selection of highly variable genes and projection of the pediatrics dataset cells into the Science cell types annotation. Be aware that this process takes large amounts of RAM memory.

```{bash}
python SLE/process_Science_pediatrics.py
```

Prepare the datasets for singleDeep.

``` bash
Rscript singleDeep/PrepareData.R --inputPath SLE/data/Science/science_processed.h5ad --fileType scanpy --sampleColumn ind_cov --clusterColumn cg_cov --clinicalColumns Status --targetColumn Status --maxCells 30000 --filterGenes --outPath SLE/data/SLE_Science

Rscript singleDeep/PrepareData.R --inputPath SLE/data/pediatrics/pediatrics_processed.h5ad --fileType scanpy --sampleColumn ind_cov --clusterColumn cg_cov --clinicalColumns Status --targetColumn Status --filterGenes --outPath SLE/data/SLE_pediatrics
```

Run the main analysis with singleDeep.

``` bash
python singleDeep/singleDeep.py --inPath SLE/data/SLE_Science --sampleColumn ind_cov --logPath SLE/log_SLE --resultsPath SLE/results_SLE/ --varColumn Status --num_epochs 500 --resultsFilenames Status --KOuter 5 --KInner 4 --saveModel
```

Predict the labels of the pediatrics dataset.

``` bash
python singleDeep/singleDeep_pretrained.py --inPath SLE/data/SLE_pediatrics --sampleColumn ind_cov --modelFile SLE/results_SLE/Status_models.pt --outFile SLE/results_SLE/validation_predictions.tsv
```

Machine learning models are used to perform the same classifications on pseudobulk data. The following command generates all the results into the corresponding folders.

```{bash}
bash SLE/run_ML_SLE.sh
```

Run the script Figure2.R to generate the Figure 2 from the article.

``` bash
RScript SLE/Figure2.R
```

Figure 2 panels will be saved in the figures folder.

## Factors influencing cell types performance

#### In silico analysis

Run the script simulation/simulation.R to generate 10 synthetic scRNA-Seq datasets, each one with 10 groups (aka cell types) with decreasing differential expression between conditions.

```{bash}
Rscript cell_types_performance/simulation.R
```

The folder simulation/data will contain the simulated datasets. The folder simulation/results_simulation_n_x will contain the singleDeep results for dataset x. Figure 3a will be saved in the figures folder.

#### Ablation analysis

Run the following script to generate subsets of the SLE Science dataset with different number of cells:

```{bash}
bash cell_types_performance/prepare_ablation.sh
```

Run singleDeep for each dataset (you may use a computational cluster to parallel this process):

```{bash}
bash cell_types_performance/run_ablation.sh
```

The results for each dataset will be stored in the ablation folder.

To generate the Figure 3 from this section results, run:

```{bash}
Rscript cell_types_performance/Figure3.R
```

Figure 3a, 3b and 3c will be saved in the figures folder.

## Dementia

Download the MTG h5ad file from <https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30> (direct link to file: <https://datasets.cellxgene.cziscience.com/77dab54a-f2a8-42fc-8c1b-3fda90622ac7.h5ad>). Rename it to *SEA_AD_MTG.h5ad* and store it into the *Dementia* folder. Now run Dementia/process_Dementia.py to select the highly variable genes and scale the data.

```{bash}
python Dementia/process_Dementia.py
```

This script will create the processed file *Dementia/Dementia_processed.h5ad*. The next step is to prepare the dataset for singleDeep.

``` bash
Rscript singleDeep/PrepareData.R --inputPath ./Dementia/Dementia_processed.h5ad --fileType scanpy --sampleColumn donor_id --clusterColumn cell_type --clinicalColumns 'disease,sex,self_reported_ethnicity,Age at death' --targetColumn disease --minCells 1000 --maxCells 30000 --filterGenes --outPath ./Dementia/data
```

Run the analyses with singleDeep for dementia prediction.

``` bash
python singleDeep/singleDeep.py --inPath ./Dementia/data/ --sampleColumn sampleID --logPath ./Dementia/log_Dementia --resultsPath ./Dementia/results_Dementia/ --varColumn disease --targetClass 0 --num_epochs 500 --resultsFilenames Status --KOuter 3 --KInner 3
```

Prepare the data for plotting:

``` bash
Rscript Dementia/preparePlot.R
```

Run the script *Dementia/Figure4.py* to generate the Figure 4 into the figures folder:

``` bash
python Dementia/Figure4.py
```

## COVID-19 (Supplementary use case)

Download the file *COVID19_ALL.h5ad.tar.gz* from <https://explore.data.humancellatlas.org/projects/5f607e50-ba22-4598-b1e9-f3d9d7a35dcc/project-matrices>. Decompress and store it into the *COVID* folder. Now run COVID/process_COVID.py to select the highly variable genes and scale the data.

```{bash}
python COVID/process_COVID.py
```

This script will create the processed file *COVID/COVID19_processed.h5ad*. The next step is to prepare the dataset for singleDeep.

``` bash
Rscript singleDeep/PrepareData.R --inputPath ./COVID/COVID19_processed.h5ad --fileType scanpy --sampleColumn sampleID --clusterColumn celltype --clinicalColumns 'CoVID-19 severity,Sample type,SARS-CoV-2,Outcome' --targetColumn 'CoVID-19 severity' --minCells 1000 --maxCells 30000 --filterGenes --outPath ./COVID/data
```

Run the analyses with singleDeep for severity and status prediction.

``` bash
python singleDeep/singleDeep.py --inPath ./COVID/data/ --sampleColumn sampleID --logPath ./COVID/log_COVID --resultsPath ./COVID/results_COVID/ --varColumn CoVID-19_severity --targetClass 2 --num_epochs 500 --resultsFilenames Severity --KOuter 3 --KInner 3
```

``` bash
python singleDeep/singleDeep.py --inPath ./COVID/data/ --sampleColumn sampleID --logPath ./COVID/log_COVIDStatus --resultsPath ./COVID/results_COVIDStatus/ --varColumn SARS-CoV-2 --targetClass 1 --num_epochs 500 --resultsFilenames Status --KOuter 3 --KInner 3
```

Run the script *COVID/SupplementaryFigure1.R* to generate the Supplementary Figure 1a and 1b into the figures folder:

``` bash
Rscript COVID/SupplementaryFigure1.R
```

## 
