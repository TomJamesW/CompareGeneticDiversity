# CompareGeneticDiversity: a framework for comparing genetic diversity measures
An R project showcasing how I compared genetic diversity measures during my MSc.

## Introduction
This project developed an experimental framework in which genetic diversity and differentiation measures can be evaluated by the proportion of time they increase when divergence or differentiation of subcommunities increases. It compares popular measures from the literature, namely $\pi$, G<sub>ST</sub>, G'<sub>ST</sub> and Jost's D,  against eachother and against a new measure I introduced, $\bar{B}$ that has now been incorporated into the package `rdiversity`.

Simulations of genetic data were carried out in the evolutionary simulation framework SLiM v.3.3, scripts for which can be found in the folder `SLiMscripts/`. A virus like genome of 8000 bases with a mutation rate of $1.2\times10^{-4}$ per base per individual per generation was used.

`deterministic.r` deterministically increases differentiation of subcommunities. Jost (Jost, 2008) showed that G<sub>ST</sub> was unfit for purpose specifically as a differentiation measure as it didn’t increase monotonically with deterministically increasing differentiation. We repeated this experiment, whilst also calculating G'<sub>ST</sub>, and $\bar{B}$ for each level of differentiation to see how it behaved. We only considered a single locus which could have many alleles. We started with two identical subcommunities, each with four equally common alleles. The original example had 1000 individuals per allele per subcommunity, but the measures are calculated using allele proportions so one individual per allele per subcommunity creates the same effect. Next, we successively added unique alleles to each subcommunity (one individual per allele) and calculated differentiation values using each of the measures. It is clear from the final plot that G<sub>ST</sub> fails this test while the rest succeed.

