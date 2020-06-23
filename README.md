# CompareGeneticDiversity: a framework for comparing genetic diversity measures
An R project showcasing how I compared genetic diversity measures during my MSc.

## Introduction
This project developed an experimental framework in which genetic diversity and differentiation measures can be evaluated by the proportion of time they increase when divergence or differentiation of subcommunities increases. It compares popular measures from the literature, namely $\pi$, F<sub>ST</sub>, G<sub>ST</sub> and G'<sub>ST</sub>,  against eachother and against a new measure I introduced, $\bar{B}$ that has now been incorporated into the package `rdiversity`.

Simulations of genetic data were carried out in the evolutionary simulation framework SLiM v.3.3, scripts for which can be found in the folder `SLiMscripts/`. A virus like genome of 8000 bases with a mutation rate of $1.2\times10^{-4}$ per base per individual per generation was used.

`deterministic.r` deterministically increases differentiation of subcommunities.