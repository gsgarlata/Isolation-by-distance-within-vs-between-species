
# Computing window-based genetic distances from a multi-species genomic dataset.


* [Compute log-nrmse between within-species and between-species isolation-by-dostance](createLOGNrmse.R): This script is used to compute ‘isolation-by-distance’ (IBD) for each genomic window by fitting a linear regression between individual pairwise genetic distances and geographic distances (logarithm). For istance, given two candidate sister species, the IBD is computed for each candidate species (within-species IBD), separately, and then by considering only individual pairwise comparisons between candidate species (between-species IBD). The fitted parameters are then used to estimate the normalized root mean squared error (NRMSE) for each genomic window (see [Van Elst et al., 2025](https://www.nature.com/articles/s41559-024-02547-w) for more details).

* [Auxiliary R functions](functions_LOGnrmse.R): This script contains all the relevant R functions needed for computing the within-species and between-species IBD as well as the NRMSE distributions.

