# Computing window-based genetic distances from a multi-species genomic dataset.

This directory provides the scripts to compute window-based genetic distances from a mulit-species genomic dataset (bcf format) based on an arbitrary window-size. It allows to select the most parsimonious window-size by comparing signal and error from the position of each individual in a PCA space across genomic windows. The signal and error variance is computed for each tested window size. This approach follows the methodology developped in [Li and Ralph, 2019](https://academic.oup.com/genetics/article/211/1/289/5931130?login=false).




