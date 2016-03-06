# MUVR
Multivariate methods with Unbiased Variable selection in R

## Issues
- fitnessOuter -> rank(colMeans) vs rank(-colMeans)
  - same min / mid / max index calculations!
- Read reg/DA from Y
- misClass -> MISS
- Recode DA
  - Fix yPred arrays for multiple Y (from RF-OlleT)
  - Fix PLS-DA using mixOmics::plsda -> plsdaInner & plsdaOuter
  - AUROC only possible for nClass=2 else MISS
- pls mode???
- Include outMod in modelReturn
- Fix vectSamp(vect,n), where length(vect)<n

## Fixed
- Validation plot (metric vs nVar)
- code rfInner.R
- code rfModOut
- check combinations of DA vs fitness (e.g. reg -> "RMSEP")

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

Buck et al 2016

Filzmoser et al 2009

Hanhineva et al 2015

Westerhuis et al 2008

## Version history
version | date | comment
:------ | :--- | :------
0.0.912 | 2016-03-06 | Added plotVAL
0.0.911 | 2016-03-03 | RF regression works in sequential mode 
0.0.910 | 2016-03-01 | Fixed PLS: DA, reg, ML, nearZeroVar
0.0.900 | 2016-02-10 | Started to work on wrapper for MV-methods. Switched from GitHub to GitLab
