# Estimating rock composition from replicate geochemical analyses: theory and application

## Description

Supporting code &amp; data for the manuscript "Estimating rock composition from replicate geochemical analyses: theory and application to magmatic rocks of the GeoPT database" by Maxime Keutgen De Greef, Gert Jan Weltje (KU Leuven) and Irène Gijbels (KU Leuven) (https://link.springer.com/article/10.1007/s11004-024-10138-5).

The main idea behind the paper is that the proper sample space for chemical analyses is the unit-hypercube, defined as the set of row-vectors whose coordinates belong to the (0,1) open set. The implication of this observation is that to do a statistical analysis of random vectors representing chemical analyses, one should first transform those random vectors living in the hypercube into their coordinates in the real space with the binary logratio (blr) function, understood as a multivariate function applying to the component of its argument the logit function :

$$ logit(x) = log(x/(1-x)) $$ 

Therefore the binary logratio (blr) function for a vector representing a chemical analyse of $D$ elements is as follows :

$$ blr(x_1,\dots,x_D) = [logit(x_1),\dots,logit(x_D)] $$

And it maps a vector of chemical analysis living in $\mathbb{H}^D$ to the $D$ dimensional real space $\mathbb{R}^D$.


## How to replicate our work 






## How to implement our methods in your own results

The optimized blr mean estimator is implemented in the first algorithm file (algo1_data_acquisition_process.R). The "optimizing.blr.mean"" is the function which computes the optimized blr mean estimate of a set of chemical analyses.

We suggest to use this function to estimate the mean of a set of chemical analyses constrained to the hypercube.

## Structure of the repository

The repository is structured as follows. 

1. The /data/ folder contains both observed chemical analyses from the GeoPT database and simulated datasets.
2. The /src/ folder contains supporting code for the paper.

Subfolders of the /data/ folder are :

1. /data/cleaned_datasets_geopt 
2. /data/datasets_with_missingvalues
3. /data/raw_datasets_geopt
4. /data/simulated_datasets


