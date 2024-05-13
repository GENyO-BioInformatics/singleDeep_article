NCellsList=("5000" "10000" "15000" "20000" "25000" "30000")
iterationList=({1..10})

for iteration in ${iterationList[@]}; do
	echo $iteration
	for NCells in ${NCellsList[@]}; do
		echo $NCells
        python singleDeep/singleDeep.py --inPath "cell_types_performance/ablation/SLE_Science_N"$NCells"_"$iteration"/" --sampleColumn ind_cov --logPath ablation/log_ablation --resultsPath "cell_types_performance/ablation/results_ablat_N"$NCells"_"$iteration"/" --varColumn Status --num_epochs 500 --resultsFilenames Status --KOuter 5 --KInner 4 --saveModel
    done
done