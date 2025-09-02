# Isolation-by-distance-within-vs-between-species

This repository provides the scripts to test whether genetic distances between candidate individuals could be explained by a model of intraspecific geographic structure. To do so, we developed a heuristic approach based on isolation-by-distance (IBD), as described in the [van Elst et al. 2025](https://www.nature.com/articles/s41559-024-02547-w) study entitled "Integrative taxonomy clarifies the evolution of a cryptic primate clade". The repository contains two folders with scripts i) for computing window-based genetic distances and ii) for computing IBD within and between candidate species. In particular:

* [Window-Based Genetic Distances](Window-Based%20Genetic%20Distances): It contains scripts to compute window-based genetic distances from a mulit-species genomic dataset based on an arbitrary window-size. It also contain scripts for selecting the most parsimonious window-size based on the lowest ratio between signal and error variance (following [Li and Ralph, 2019](https://academic.oup.com/genetics/article/211/1/289/5931130?login=false) approach).

* [Window-Based IBD within-vs-between species](Window-Based%20IBD%20within-vs-between%20species): It contains scripts to compute isolation-by-distance within and between candidates by correlating individual genetic distances with geographic distances (on log scale) and quantify deviations of observed genetic distances between candidates from those predicted by the within-candidate geographic clines in genetic distance.

For any question and doubts send an email to gabriele.sgarlata[at]gmail[dot]com.
