for i in list(range(1,11)):
  for j in ['model','model-fam','model-fam-geo']:
    f = open('invoker_{}_{}.sh'.format(j,i),'w')
    print('Rscript run_model.R {} {}'.format(j,i),file=f)
    f.close()
    f = open('run_{}_{}.sh'.format(j,i),'w')
    print("""#!/bin/bash
#SBATCH --qos=medium
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

module load generic
module load singularity

srun singularity exec -u ~/r_latest.sif bash invoker_{}_{}.sh""".format(j,i),file=f)
    f.close()