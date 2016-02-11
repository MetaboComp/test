# MUVR
Multivariate methods with Unbiased Variable selection in R

## Issues
- fix vectSamp(vect,n), where length(vect)<n
- code rfInner.R
- code rfModOut
- Fix yPred arrays for multiple Y (from RF-OlleT)
- check combinations of DA vs fitness (e.g. reg -> "RMSEP")
- Include modReturn in modelReturn

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
0.0.900 | 2016-02-10 | Started to work on wrapper for MV-methods. Switched from GitHub to GitLab
