
module load new
module load /cluster/apps/imsb/modules
module load r/3.5.1
module load mpfr
module load open_mpi

cd ~/mysonas/PRPF8/analysis

bsub -J genomic_coord -R "rusage[mem=5000,scratch=5000]" -W 48:00 Rscript --vanilla ./CCprofilerAnalysis/final/02_tracesProcessing.R
