
module load new
module load /cluster/apps/imsb/modules
module load r/3.5.1
module load mpfr
module load open_mpi

cd ~/mysonas/PRPF8/analysis

# bsub -J 01_importData -R "rusage[mem=20000,scratch=20000]" Rscript --vanilla ./CCprofilerAnalysis/update/01_importData.R
bsub -J 02_tracesProcessing -R "rusage[mem=20000,scratch=20000]" -W 24:00 Rscript --vanilla ./CCprofilerAnalysis/update/02_tracesProcessing.R
bsub -J 03_annotateProteoforms -w "02_tracesProcessing" -R "rusage[mem=20000,scratch=20000]" Rscript --vanilla ./CCprofilerAnalysis/update/03_annotateProteoforms.R

bsub -R "rusage[mem=20000,scratch=20000]" -W 24:00 Rscript --vanilla ./CCprofilerAnalysis/benchmark/05_proteinFeatures_proteoformResolved_update.R
bsub -R "rusage[mem=20000,scratch=20000]" -W 4:00 Rscript --vanilla ./CCprofilerAnalysis/benchmark/06_proteinFeatures_proteoformResolved_resolved.R
