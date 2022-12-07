#' Calculate permutation p-value of actual model performance vs null hypothesis distribution
#'
#' `pPerm` will calculate the cumulative (1-tailed) probability of `actual` belonging to `permutation_distribution`.
#' Side is guessed by actual value compared to median(permutation_distribution).
#' Test is performed on original data OR ranked for non-parametric statistics.
#' @param actual Actual model performance (e.g. misclassifications or Q2)
#' @param permutation_distribution Null hypothesis distribution from permutation test (same metric as `actual`)
#' @param side Smaller or greater than (automatically guessed if omitted) (Q2 and AUC is a "greater than" test, whereas misclassifications is "smaller than")
#' @param type Standard Student's t distribution ('t') or Student's t on rank-transformed data for nonparametric test ('non')
#' @return p-value
#' @export
pPerm=function(actual,                             ###a value
               permutation_distribution,            ###a distribution
               side=c('smaller','greater'),
               type=c('t','non',"smooth","ecdf","rank")) {
  p_object<-list()
##########################################################################################################################
#it needs to be take into consideration when the type of side is error
  if(is.numeric(actual)==F){
    stop("actual needs to be numeric")
  }
  if(is.numeric(permutation_distribution)==F){
    stop("permutation_distribution needs to be a numeric distribution")
  }
  if(length(permutation_distribution)<5){
    stop("permutation_distribution has too few values to form a distribution")
  }
  if(!missing(type)){if(!type%in%c('t','non',"smooth","ecdf","rank")){
    stop("This type can not be implemented")
  }
    }
  if(missing(type)) {type='t'}     ###Student's t distribution
  if(!missing(side)){if(side!="smaller"&side!="greater"){
    stop("This side can not be implemented")
  }
    }
  if(missing(side)) {side=ifelse(actual<median(permutation_distribution),
                                 'smaller',
                                 'greater')}   ###This is to see which side

######################################################################################################################################################################################################

  if (type=='ecdf') {
    if(side=="smaller"){
      sort_x<-sort(permutation_distribution)
      fun_ecdf<-ecdf(sort_x)
      ecdf_curve <- fun_ecdf(sort_x)

      if(actual<min(sort_x)){
        p_actual<-"<0.01"
      }else{
        for(i in 1:(length(permutation_distribution)-1)){
          if(actual>sort_x[i]&actual<=sort_x[i+1]){
            p_actual<-ecdf_curve[i]
          }

        }
        if(p_actual==0){p_actual<-"<0.01"}
      }
      p<-p_actual

   }else{

     sort_x<-sort(permutation_distribution)
     fun_ecdf<-ecdf(sort_x)
     ecdf_curve <- fun_ecdf(sort_x)

     if(actual>max(sort_x)){
       p_actual<-"<0.01"
     }else{
       for(i in 1:(length(permutation_distribution)-1)){
         if(actual>=sort_x[i]&actual<sort_x[i+1]){
           p_actual<-1-ecdf_curve[i]
         }
       }
       if(p_actual==0){p_actual<-"<0.01"}
     }
     p<-p_actual

   }
    p_object$p<-p
    p_object$points<-NULL
  }



  ######################################################################################################################
  ## redo p value from cumulative curve
  if (type=='non') {
    if(side=="smaller"){
     rank=rank(c(actual,permutation_distribution))     ###the sequence of each value
     actual=rank[1]              ###actual is not the smallest if side is greater this apply to Q2 and AUC
     permutation_distribution=rank[-1]
    }else{rank=rank(c(actual,permutation_distribution))
          actual=rank[(length(permutation_distribution)+1)]
          permutation_distribution=rank[-(length(permutation_distribution)+1)]
    }

    p=pt((actual-mean(permutation_distribution))/sd(permutation_distribution),    ##### pt gives the probability before the input point

         df=length(permutation_distribution)-1)

    if (side=='greater') {p=1-p }
    p_object$p<-p
    p_object$points<-NULL
  }

  ######################################################################################################################
  if (type=='rank') {
    if(side=="smaller"){
      rank=rank(c(actual,permutation_distribution))     ###the sequence of each value
      actual=rank[1]/length(rank)              ###actual is not the smallest if side is greater this apply to Q2 and AUC
      permutation_distribution=rank[-1]
    }else{rank=rank(c(actual,permutation_distribution))
    actual=1/length(rank)+1-(rank[1]/length(rank))
    permutation_distribution=rank[-1]
    }

    p<-actual
    cat("/n Resolution limited to ",1/length(rank))
    p_object$p<-p
    p_object$points<-NULL

  }

###############################################################################################################################################################################################

  if(type=="t"){
  p=pt((actual-mean(permutation_distribution))/sd(permutation_distribution),    ##### pt gives the probability before the input point

       df=length(permutation_distribution)-1)

  if (side=='greater') {p=1-p }    ####this is to solve which side p value is

  p_object$p<-p
  p_object$points<-NULL
  }


  if(type=="smooth"){
    e = 0.01 * diff(range(permutation_distribution))
    if(actual>=max(permutation_distribution)){

    from=min(permutation_distribution)-(actual-max(permutation_distribution))-e
    to=actual+e
    }
    if(actual<=min(permutation_distribution)){
      to=min(permutation_distribution)-actual+max(permutation_distribution)+e
      from=actual-e
    }

    dens = density(permutation_distribution,
                   adjust=0.01,
                   from=from,
                   to=to,
                   n = 10000)

    if(side=="smaller"){
      sort_x<-sort(dens$x)
      fun_ecdf<-ecdf(sort_x)
      ecdf_curve <- fun_ecdf(sort_x)


        for(i in 1:(length(sort_x)-1)){
          if(actual>sort_x[i]&actual<=sort_x[i+1]){

            p_actual<-ecdf_curve[i]
          }
        }


      p<-p_actual

    }else{

      sort_x<-sort(dens$x)
      fun_ecdf<-ecdf(sort_x)
      ecdf_curve <- fun_ecdf(sort_x)

        for(i in 1:(length(sort_x)-1)){
          if(actual>=sort_x[i]&actual<sort_x[i+1]){
            p_actual<-1-ecdf_curve[i]
          }
        }


      p<-p_actual

    }

  p_object$p<-p
  p_object$points<-cbind.data.frame(sort_x,ecdf_curve)
  p_object$dens<-dens
  }


  return(p_object)
}
