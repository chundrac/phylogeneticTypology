python3 config_slurm.py

module use /sapps/etc/modules/start/
module load generic

for j in 'model' 'model-fam' 'model-geo-fam' 'model-geo'
do
	for i in {1..10}
	do     
	export SINGULARITY_BIND="/home/cluster,/net/cephfs/home,/data,/net/cephfs/data,/scratch,/net/cephfs/scratch,/net/cephfs/shares,/sapps,/net/cephfs/sapps"
	sbatch run\_$j\_$i.sh
done
done

rm *.sh