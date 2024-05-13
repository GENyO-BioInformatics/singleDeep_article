NCellsList=("5000" "10000" "15000" "20000" "25000" "30000")
iterationList=({1..10})

for iteration in ${iterationList[@]}; do
	echo $iteration
	for NCells in ${NCellsList[@]}; do
		echo $NCells
		Rscript singleDeep/PrepareData.R --inputPath SLE/data/Science/science_processed.h5ad --fileType scanpy --sampleColumn ind_cov --clusterColumn cg_cov --clinicalColumns Status --targetColumn Status --minCells 25000 --maxCells $NCells --filterGenes --outPath "cell_types_performance/ablation/SLE_Science_N"$NCells"_"$iteration --seed $iteration
	done
done
