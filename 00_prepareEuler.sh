
module load new
module load /cluster/apps/imsb/modules
module load r/3.4.0
module load mpfr
module load open_mpi

cd ~/mysonas/PRPF8/analysis

# do script one manually because of biomart issues
# bsub -J 01_importData Rscript --vanilla ./01_importData.R
# bsub -J 02_tracesProcessing -w "01_importData" Rscript --vanilla ./02_tracesProcessing.R

bsub -J 02_tracesProcessing Rscript --vanilla ./02_tracesProcessing.R
bsub -J 03_tracesStats -w "02_tracesProcessing" Rscript --vanilla ./03_tracesStats.R
bsub -J 04_plotExamplePepTraces -w "03_tracesStats" Rscript --vanilla ./04_plotExamplePepTraces.R
bsub -J 05_proteinQuant -w "04_plotExamplePepTraces" Rscript --vanilla ./05_proteinQuant.R
bsub -J 06_plotExampleProteins -w "05_proteinQuant" Rscript --vanilla ./06_plotExampleProteins.R

# run DAVID
bsub -J 07_massDistributionDavid Rscript --vanilla ./07_massDistributionDavid.R

# protein feature finding
# bsub -n 10 -J 08_proteinFeatureFinding_plus Rscript --vanilla ./08_proteinFeatureFinding_plus.R
# bsub -n 10 -J 08_proteinFeatureFinding_minus Rscript --vanilla ./08_proteinFeatureFinding_minus.R
# bsub -n 10 -J 08_proteinFeatureFinding Rscript --vanilla ./08_proteinFeatureFinding_all.R
# bsub -J 10_proteinFeatureStats Rscript --vanilla ./10_proteinFeatureStats.R
# integrated protein feature finding
bsub -n 10 -J 08_proteinFeatureFinding_all Rscript --vanilla ./CCprofilerAnalysis/08_proteinFeatureFinding_all.R

# complex feature finding
bsub -n 10 -J 09_complexFeatureFinding Rscript --vanilla ./09_complexFeatureFinding.R
bsub -J 11_complexFeatureStats -w "09_complexFeatureFinding" Rscript --vanilla ./11_complexFeatureStats.R


bsub -J 05_proteinQuant -R "rusage[mem=500000,scratch=500000]" -W 4:00 Rscript --vanilla ./CCprofilerAnalysis/05_proteinQuant.R
