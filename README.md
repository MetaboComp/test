# MUVR
Multivariate methods with Unbiased Variable selection in R

## General description
This is a package dedicated to predictive multivariate modelling for metabolomics.
- Types: Classification and regression
- Model cores: PLS and Random Forest
- Data structures: Paired and unpaired
- Validation: rdCV (Westerhuis et al 2008, Filzmoser et al 2009)
- Variable selection: Performed internally.  
  The unbiased VS stems from being tuned in the inner CV loop.  

## Installation
Install `devtools` to be able to install packages from Git repositories.

Install `MUVR` package by:

`devtools::install_git("https://gitlab.com/CarlBrunius/MUVR.git")`

In addition to functions relevant for within/between batch correction, data is provided to accurately reproduce figures from the original *Brunius et al* paper (see below).

## Examples of use
MV-core       | response       | data structure | comment
:------       | :-------       | :------------- | :------
PLS           | regression     | unpaired       | Used in Hanhineva et al 2015 on LC-MS metabolomics data (linear code; before MUVR)
PLS           | multilevel     | paired         | Used in 'BioDiVa', 'Satiety' and 'CRISP' projects.
Random Forest | classification | unpaired       | Used in Buck et al 2016 on 16S microbiomics data (linear code; before MUVR)
Random Forest | multilevel     | paired         | Used in 'BioDiVa' project.

## References

- *Buck M, Nilsson LKJ, Brunius C, Dabire RK, Hopkins R, Terenius O, 2016. Bacterial associations reveal spatial population dynamics in Anopheles gambiae mosquitoes. Scientific Reports; Accepted.*
- *Filzmoser P, Liebmann B, Varmuza K, 2009. Repeated double cross validation. Journal of Chemometrics 23(4), 160-171.*
- *Hanhineva K, Brunius C, Andersson A, Marklund M, Juvonen R, Keski-Rahkonen P, Auriola S, Landberg R., 2015. Discovery of urinary biomarkers of whole grain rye intake in free-living subjects using nontargeted LC-MS metabolite profiling, Molecular Nutrition and Food Research 59(11), 2315-25.*
- *Westerhuis JA, Hoefsloot HCJ, Smit S, Vis DJ, Smilde AK, Velzen EJJ, Duijnhoven JPM, Dorsten FA, 2008. Assessment of PLSDA cross validation. Metabolomics 4(1), 81-89.*

## Version history
version | date | comment
:------ | :--- | :------
0.0.965 | 2017-11-06 | Segmentation for DA & dependent samples -> vectSamp(unikID)
0.0.964 | 2017-10-16 | Fixed geometric averaging for nVar in MUVR, fixed plotStability, plotVAL and workflows in inst/workflow folder
0.0.963 | 2017-10-16 | Fixes in plotStability, plotVAL and MUVR2
0.0.962 | 2017-10-04 | Fixed misindexing in pls function (nzv in now before variable allocation)
0.0.961 | 2017-10-03 | Removed dependency on mixOmics by adapting pls functions from mixOmics v 5.2.0. Revamped MUVR.
0.0.960 | 2017-09-22 | Fixed bug with PLS-DA and no successful model -> return NA. Added `cut` to plotVIP and minor fix in MUVR
0.0.959 | 2017-06-20 | Fixed plsInner error with DA and matrix format
0.0.958 | 2017-06-01 | Fixed plsInner 'drop' error
0.0.957 | 2017-05-22 | Added `getVIP`
0.0.956 | 2017-05-18 | removed old nzv code in plsInner and MUVR
0.0.955 | 2017-04-20 | `plsInner` now decreases ncomp if producing NAs. Also cleaned up some code in `MUVR`
0.0.954 | 2017-03-09 | Added `biplotPLS`
0.0.953 | 2017-02-14 | Added `pPerm`, `plotPerm` and `plotVIP`
0.0.952 | 2017-02-07 | Updated plsInner -> removes comp if error
0.0.951 | 2016-12-12 | Updated plsInner -> deals with NA
0.0.950 | 2016-11-16 | Added first version of predict function `predMV` (for PLS regression) and tweaked MUVR to allow for submodel predictions
0.0.945 | 2016-11-14 | Added choice of Y in ML (for easier permutation analysis)
0.0.944 | 2016-11-10 | Added choice (non)-parallel
0.0.943 | 2016-07-25 | MVWrap -> MUVR & Included VIPPerRep in MUVR and RDCV.
0.0.940 | 2016-06-08 | Added `rdCV` algorithm. Updated `workflow.r`
0.0.932 | 2016-04-27 | Updated `rfInner` and `MVWrap` for ML-RF.
0.0.931 | 2016-03-31 | Added R2 and Q2 for regression. Fixed yPred (previously Min for all methods!)
0.0.930 | 2016-03-08 | Recoded ML to regression. Added plotMV. Migrated `Issues` to GitLab.
0.0.920 | 2016-03-07 | Recoded DA (also multiclass) for PLS and RF
0.0.913 | 2016-03-06 | New index functions: apply/rank -> rank(colMeans)
0.0.912 | 2016-03-06 | Added plotVAL
0.0.911 | 2016-03-03 | RF regression works in sequential mode 
0.0.910 | 2016-03-01 | Fixed PLS: DA, reg, ML, nearZeroVar
0.0.900 | 2016-02-10 | Started to work on wrapper for MV-methods. Switched from GitHub to GitLab
