# Computing window-based genetic distances from a multi-species genomic dataset.

This directory provides the scripts to compute window-based genetic distances from a mulit-species genomic dataset (bcf format) based on an arbitrary window-size. It allows to select the most parsimonious window-size by comparing signal and error from the position of each individual in a PCA space across genomic windows. The signal and error variance is computed for each tested window size. This approach follows the methodology developped in [Li and Ralph, 2019](https://academic.oup.com/genetics/article/211/1/289/5931130?login=false).

*[Compute PCA for each genomic window](createPCs_list_by_chromosome.R): This script is used to create a file for each chromosome, which containes the PCA positioning of each individual across the various genomic windows contained in the chromosome (given a defined window-size). For the purposes of selecting the most-parimonious window-size (higher signal-to-error variance), one should run this script for various window-sizes.







