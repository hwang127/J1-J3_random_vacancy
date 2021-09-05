#!/bin/csh


foreach nsamples (20000)
foreach seed ( ` seq 1 10 ` )
foreach nwait(0)
foreach L (76 )# 76 95 114 133 )
foreach xi (100 )  #these are actually the num
foreach T (0.399 0.4 0.4005 0.401 0.402 0.403)#(0.4 0.401 0.402 0.4025 0.403 0.404 0.405)#(0.42 0.425 0.426 0.427 0.428 0.429 0.43)#(0.42 0.415 0.41 0.407 0.405 0.404 0.402 0.4 0.395 0.39)#0.385 0.38 0.375 0.37 0.365 0.36 0.355 0.35)
foreach alpha(0)
foreach beta(0)#9 36 100  will be taken sqrta
foreach outNum (1)
foreach J(-0.6165)
foreach Jthird(3.0)

# Input file
cat >/home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/jobs/in_seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}.txt << EOFm
nsamples=${nsamples}
nwait=0
L=${L}
J=${J}
Jthird=${Jthird}
T=${T}
beta=${beta}
seed=${seed}
outNum=${outNum}
outfilename=/home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/runs/seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}_

EOFm

cat >/home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/jobs/job_seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}<<EOFm
#!/bin/bash -l
# Batch Queue Script
#SBATCH --time=72:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=hwang127@jhu.edu
#SBATCH --cpus-per-task=1
#SBATCH --account=olegt
/home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/mcfile /home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/jobs/in_seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}.txt /home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/1_${L}_${xi}_sample.txt 

EOFm

chmod 755 /home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/jobs/job_seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}
sbatch /home-2/hwang127@jhu.edu/scratch/newUpdate/vacancy/energy/jobs/job_seed_${seed}_T_${T}_beta_${beta}_xi_${xi}_L_${L}
end
end
end
end
end
end
end
end
end
end
