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

algorithm | MV-core | response | data structure | comment
:-------- | :------- | :------------- | :------ | :------
PLS-VS | PLS | class/reg | unpaired | Parallellised. Used in Hanhineva et al 2015.
ML-PLS-VS | PLS | class/reg | paired (multilevel) | Parallellised. Used in 'BioDiVa' and 'Satiety' projects.
RF-VS | Random Forest | classification | unpaired | Parallellised. Used in Buck et al 2016.
ML-RF-VS | Random Forest | classification | paired (multilevel) | Parallellised. Used in 'BioDiVa' project.

## References

- *Buck M, Nilsson LKJ, Brunius C, Dabire RK, Hopkins R, Terenius O, 2016. Bacterial associations reveal spatial population dynamics in Anopheles gambiae mosquitoes. Scientific Reports; Accepted.*
- *Filzmoser P, Liebmann B, Varmuza K, 2009. Repeated double cross validation. Journal of Chemometrics 23(4), 160-171.*
- *Hanhineva K, Brunius C, Andersson A, Marklund M, Juvonen R, Keski-Rahkonen P, Auriola S, Landberg R., 2015. Discovery of urinary biomarkers of whole grain rye intake in free-living subjects using nontargeted LC-MS metabolite profiling, Molecular Nutrition and Food Research 59(11), 2315-25.*
- *Westerhuis JA, Hoefsloot HCJ, Smit S, Vis DJ, Smilde AK, Velzen EJJ, Duijnhoven JPM, Dorsten FA, 2008. Assessment of PLSDA cross validation. Metabolomics 4(1), 81-89.*

## Version history
version | date | comment
:------ | :--- | :------
0.0.930 | 2016-03-08 | Recoded ML to regression. Added plotMV. Migrated `Issues` to GitLab.
0.0.920 | 2016-03-07 | Recoded DA (also multiclass) for PLS and RF
0.0.913 | 2016-03-06 | New index functions: apply/rank -> rank(colMeans)
0.0.912 | 2016-03-06 | Added plotVAL
0.0.911 | 2016-03-03 | RF regression works in sequential mode 
0.0.910 | 2016-03-01 | Fixed PLS: DA, reg, ML, nearZeroVar
0.0.900 | 2016-02-10 | Started to work on wrapper for MV-methods. Switched from GitHub to GitLab
