# Computing window-based genetic distances from a multi-species genomic dataset.

This directory provides the scripts to compute window-based genetic distances from a mulit-species genomic dataset (bcf format) based on an arbitrary window-size. It allows to select the most parsimonious window-size by comparing signal and error from the position of each individual in a PCA space across genomic windows. The signal and error variance is computed for each tested window size. This approach follows the methodology developped in [Li and Ralph, 2019](https://academic.oup.com/genetics/article/211/1/289/5931130?login=false).

* [Compute PCA for each genomic window](createPCs_list_by_chromosome.R): This script is used to create a file for each chromosome, which containes the PCA positioning of each individual across the various genomic windows contained in the chromosome (given a defined window-size). For the purposes of selecting the most-parimonious window-size (higher signal-to-error variance), one should run this script for various window-sizes.

* [Compute Signal-to-Error](computeSignal2ErrorPCs.R): This script is used to quantify meaningful variation in genetic patterns across the genome (signal) and variation that represent demographic noise (noise). ‘signal’ is computed as the standard deviation of the sample’s position on PC1 over all windows, averaged over samples. ‘error’ is computed as the standard error of the sample’s position on PC1 across windows and samples using the block jackknife. See [Li and Ralph, 2019](https://academic.oup.com/genetics/article/211/1/289/5931130?login=false) for more details. For the purposes of selecting the most-parimonious window-size (higher signal-to-error variance), one should run this script for various window-sizes.

* [Auxiliary R functions for Signal-to-Error analysis](functions_signal2error_PCs.R): This script contains all the relevant R functions needed for computing the standard deviation of the sample’s position on PC1 over all windows (signal) and the standard error of the sample’s position on PC1 across windows and samples using the block jackknife (noise). This script is loaded in [Compute Signal-to-Error](computeSignal2ErrorPCs.R) R script.

* [Plot Signal-to-Error across Window-Sizes](makeSignal2ErrorDELTA_PLOT.R): This script is used to plot the difference in signal and noise across window sizes.

* [Compute Window-Based Genetic Distances](createGeneticDistancesPIXY_like.R): This script is used to compute ‘Isolation-by-distance’ between candidate species. It first split the multispecies genomic dataset into windows of fixed size. Note that the window size is determined by number of SNPs and not by number of base-pairs. This means that a window size of 1000 will contain 1000 SNPs which may span over, let's say, a genomic segment of 2000 bp.

* [Auxiliary R functions for Window-Based Genetic Distances](functions_regressionPIXYlike.R): This script contains all the relevant R functions needed for computing individual pairwise genetic differences. This script is loaded in the [Compute Window-Based Genetic Distances](createGeneticDistancesPIXY_like.R) R script.








