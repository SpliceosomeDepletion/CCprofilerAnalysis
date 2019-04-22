module load new
module load /cluster/apps/imsb/modules
module load r/3.5.1
module load mpfr
module load open_mpi

cd ~/myspectrumscale/

# bsub -J genomic_coord -R "rusage[mem=5000,scratch=5000]" -W 48:00 Rscript --vanilla ./PRPF8/analysis/CCprofilerAnalysis/thesis/02B_annotateCoordinates.R 
