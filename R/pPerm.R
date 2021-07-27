#' Calculate permutation p-value of actual model performance vs null hypothesis distribution
#'
#' `pPerm` will calculate the cumulative (1-tailed) probability of `actual` belonging to `h0`.
#' Side is guessed by actual value compared to median(h0).
#' Test is performed on original data OR ranked for non-parametric statistics.
#' @param actual Actual model performance (e.g. misclassifications or Q2)
#' @param h0 Null hypothesis distribution from permutation test (same metric as `actual`)
#' @param side Smaller or greater than (automatically guessed if omitted) (Q2 and AUC is a "greater than" test, whereas misclassifications is "smaller than")
#' @param type Standard Student's t distribution ('t') or Student's t on rank-transformed data for nonparametric test ('non')
#' @return p-value
#' @export
pPerm=function(actual,
               h0,
               side=c('smaller','greater'),
               type=c('t','non')) {
##########################################################################################################################
#it needs to be take into consideration when the type od side is error
  if(is.numeric(actual)==F)stop("actual needs to be numeric")
  if(is.numeric(h0)=F)stop("h0 needs to be a numeric distribution")
  if(length(h0)<5)stop("h0 has too view values to form a distribution")
  if(!missing(type)){if(type!="t"&type!="non")stop("This type can not be implemented")}
  if(missing(type)) type='t'     ###Student's t distribution
  if(!missing(side)){if(side!="smaller"&side!="greater")stop("This side can not be implemented")}
  if(missing(side)) side=ifelse(actual<median(h0),'smaller','greater')   ###This is to see which side

######################################################################################################################################################################################################

  if (type=='non') {
    if(side=="smaller"){
     rank=rank(c(actual,h0))     ###the sequence of each value
     actual=rank[1]              ###actual is not the smallest if side is greater this apply to Q2 and AUC
     h0=rank[-1]
    }else{rank=rank(c(actual,h0))
          actual=rank[(length(h0)+1)]
          h0=rank[-(length(h0)+1)]
          }
  }
###############################################################################################################################################################################################
  p=pt((actual-mean(h0))/sd(h0),    ##### pt gives the probability before the input point

       df=length(h0)-1)

  if (side=='greater') p=1-p     ####this is to solve which side p value is

  return(p)
}
