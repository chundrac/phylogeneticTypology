for j in 'model' 'model-fam' 'model-geo-fam' 'model-geo'
do
	for i in {1..25}
	Rscript run_model $j $i
done