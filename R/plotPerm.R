#' Plot for comparison of actual model fitness vs permutation
#'
#' Plots histogram of null hypothesis (permutation) distribution, actual model fitness and cumulative p-value.
#' Plot defaults to "greater than" or "smaller than" tests and cumulative probability in Student's t-distribution
#'
#' @param actual Actual model fitness (e.g. Q2, AUROC or number of misclassifications)
#' @param distribution Null hypothesis (permutation) distribution of similar metric as `actual`
#' @param xlab Label for x-axis (e.g. 'Q2 using real value',"Q2 using distributions","BER" 'AUROC', or 'Misclassifications')
#' @param ylab label for y-axis
#' @param side Cumulative p either "greater" or "smaller" than H0 distribution (defaults to side of median(H0))
#' @param type c('t','non',"smooth","rank","ecdf")
#' @param xlim Choice of user-specified x-limits (if default is not adequate)
#' @param ylim Choice of user-specified y-limits (if default is not adequate)
#' @param breaks Choice of user-specified histogram breaks (if default is not adequate)
#' @param main Choice of user-specified plot title
#' @param permutation_visual choice of showing median or mean or none
#' @param pos Choice of position of p-value label (if default is not adequate)
#' @param curve if add curve or not base on the mid
#' @param extend how many percenrtage of the orignical range do we start
#' @param show_p if p value is added to the figure
#' @param multiple_p_shown show many p values
#' @param round_number How many digits does it keep
#' @return Plot
#' @export
plotPerm=function(actual,
                  distribution,     ####a distribution
                  xlab=NULL,
                  side=c('greater','smaller'),
                  type="t",
                  ylab=NULL,
                  #               type=c('t','non',"smooth","ecdf","rank"),
                  xlim,
                  ylim=NULL,
                  breaks='Sturges',
                  pos,  ####Choice of position of p-value label (if default is not adequate)
                  main=NULL,
                  permutation_visual="none",
                  curve=F,
                  extend=0.1,
                  multiple_p_shown=NULL,
                  show_p=T,
                  round_number=4
                  ) {

  if(!is.null(multiple_p_shown)){
    if(!any(multiple_p_shown%in%c('t','non',"smooth","ecdf","rank"))){
    stop("This type can not be implemented")
    }}
  if(!permutation_visual%in%c("mean","median","none")){stop("this type not supoorted")}
  if(!missing(type)){if(!any(type%in%c('t','non',"smooth","ecdf","rank"))){stop("This type can not be implemented")}}
  if(missing(type)) {type=='t'}

  if(!missing(side)){if(side!="smaller"&side!="greater")stop("This side can not be implemented")}

######when it is Q2 or MISS
  if(!missing(pos)){pos_rev=6-pos}
  if(missing(side)) {side=ifelse(actual<median(distribution),'smaller','greater')}

  if(missing(pos)) {pos<-ifelse(side=='smaller',4,2)
                    pos_rev<-ifelse(side=='smaller',2,4)}
  ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
  ##respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
if(!length(multiple_p_shown)>=2){
  pP=pPerm(actual,
           distribution,
           side,
           type=type,
           extend=extend)     ####calculate p value
  #if(missing(ylim)) {
  #  xlim =c(0,max(h$density))
  #}

    ran<-range(c(actual,distribution))
  from=ran[1]-diff(ran)*extend
  to=ran[2]+diff(ran)*extend
  if(missing(xlim)) {
    xlim =c(from,to)
  }


    h=hist(distribution,
          breaks,
          xlim=xlim,
          #ylim=ylim,
          axes=F,     ###remove both axes
          xlab=xlab,
          ylab=ylab,
          yaxt='n',
          freq=FALSE,   ##### if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one).
          main=main)

  h2=max(h$density)*.75  ######as estimated density values, This is to decide how high the vertical line will be drawn
  if(curve==T){
    if(type=="smooth"){
      ran<-range(c(actual,distribution))
      from=ran[1]-diff(ran)*extend
      to=ran[2]+diff(ran)*extend

    #dx=density(distribution,
    #           #adjust=0.0001,
    #           from=from,
    #           to=to,
    #           n = 100000)
    #dx_line=density(distribution,
    #           #adjust=0.0001,
    #           from=xlim[1],
    #           to=xlim[2],
    #           n = 100000)

    lines(pP$dens,lwd = 2, col = "red")
    }

    if(type=="t"){
      #e = 1 * diff(range(distribution))

      x_values<-distribution
      x_values_sim <- seq(min(distribution), max(distribution), length = 1000)
    #  curve(dt((x-mean(x_values))/(sd(x_values)/sqrt(length(x_values))),
    #               df=length(x_values)-1))
      y_values<-dt((x_values_sim-mean(x_values))/(sd(x_values)/sqrt(length(x_values))),
                   df=length(x_values)-1)
      #plot(x_values_sim, y_values*length(x_values))
      lines(x_values_sim, y_values,#*length(x_values_sim),
            lwd = 2,col = "red")

    }
  }

  axis(1,pos=0)   ###the coordinate at which the axis line is to be drawn: if not NA this overrides the value of line.

  ####comment out if you don't want y axis
 # if(side=='smaller') {axis(2,pos=0,las=1)   #### the style of axis labels. (0=parallel, 1=all horizontal, 2=all perpendicular to axis, 3=all vertical)
#  }else {axis(2,pos=h$breaks[1],las=1)}

  lines(rep(actual,2),     ###x1,x2 for the line
        c(0,h2))          ##y1 ,y2 forthe line

if(show_p==T){
  if(!is.nan(pP$p)&is.numeric(pP$p)){
  text(actual,    ###x position of the text
       h2,        ##y position of the text
       pos=pos,
       ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
       ##respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
       labels=paste('p=',
                    signif(pP$p,round_number),
                    sep='')) ####what is the text
       ##For signif the recognized values of digits are 1...22, and non-missing values are rounded to the nearest integer in that range.
       ##Complex numbers are rounded to retain the specified number of digits in the larger of the components.
       ##Each element of the vector is rounded individually, unlike printing.
  }

  if(!is.nan(pP$p)&!is.numeric(pP$p)){
    text(actual,    ###x position of the text
         h2,        ##y position of the text
         pos=pos,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         ##respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
         labels=paste("p",pP$p,sep='')) ####what is the text
    ##For signif the recognized values of digits are 1...22, and non-missing values are rounded to the nearest integer in that range.
    ##Complex numbers are rounded to retain the specified number of digits in the larger of the components.
    ##Each element of the vector is rounded individually, unlike printing.
  }

  if(!is.null(pP$points)){
    ################################ add curves

    #lines(dens)
  }

}
    text(actual,    ###x position of the text
       max(h$density)*0.003,        ##y position of the text
       pos=3,
       ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
       labels=paste0(signif(actual,round_number))
  )
if(permutation_visual=="mean"){
  text(median(distribution),    ###x position of the text
       max(h$density)*0.003,        ##y position of the text
       pos=3,
       ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
       labels=paste0("mean=",
                    signif(mean(distribution),round_number)
                    )
  )
}

  if(permutation_visual=="median"){
    text(median(distribution),    ###x position of the text
         max(h$density)*0.003,        ##y position of the text
         pos=3,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         labels=paste0("median=",
                       signif(median(distribution),round_number)
         )
    )
  }

  if(permutation_visual=="none"|missing(permutation_visual)){

  }
}else {
  #############################################################################################
  ran<-range(c(actual,distribution))
  from=ran[1]-diff(ran)*extend
  to=ran[2]+diff(ran)*extend
  if(missing(xlim)) {
    xlim =c(from,to)
  }
  h=hist(distribution,
         breaks,
         xlim=xlim,
        # ylim=ylim,
         axes=F,     ###remove both axes
         xlab=xlab,
        ylab=ylab,
         yaxt='n',
        freq=FALSE,   ##### if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one).
         main=main)

  ppp<-list()
  pP<-c()
  for(i in 1:length(multiple_p_shown)){
  ppp[[i]]=pPerm(actual,
           distribution,
           side,
           type=multiple_p_shown[i],
           extend=extend)     ####calculate p value
  pP[i]<-ppp[[i]]$p
  }



    if("t"%in%multiple_p_shown&!"smooth"%in%multiple_p_shown){
      #e = 1 * diff(range(distribution))


      x_values<-distribution
      x_values_sim <- seq(min(distribution), max(distribution), length = 1000)
      y_values<-dt((x_values_sim-mean(x_values))/(sd(x_values)/sqrt(length(x_values))),
                   df=length(x_values)-1)
      #plot(x_values_sim, y_values*length(x_values))
      lines(x_values_sim, y_values*length(x_values), lwd = 2,col = "green")
      #
      #  y_values <- dt(x_values,df=length(distribution)-1)
      #y_values <- y_values * diff(h$mids[1:2]) * length(distribution)

  #    x_values <- seq(from, to, length = 1000)
   #   y_values <- dt(x_values,df=length(distribution)-1)

      legend('topright',
             legend=c("t"),
             lty=1,
             cex =0.5,
             trace = F,####line type
             col="green",
             bty='n')

    }

  h2=max(h$density)*.75  ######as estimated density values, This is to decide how high the vertical line will be drawn

    if("smooth"%in%multiple_p_shown&!"t"%in%multiple_p_shown){
      lines(ppp[[which(multiple_p_shown=="smooth")]]$dens,lwd = 2, col = "red")
      legend('topright',
             legend=c("smooth"),
             lty=1,
             cex =0.5,
             trace = F,####line type
             col="red",
             bty='n')
    }
  if("smooth"%in%multiple_p_shown&"t"%in%multiple_p_shown){
    x_values<-distribution
    x_values_sim <- seq(min(distribution), max(distribution), length = 1000)
    y_values<-dt((x_values_sim-mean(x_values))/(sd(x_values)/sqrt(length(x_values))),
                 df=length(x_values)-1)
    #plot(x_values_sim, y_values*length(x_values))
    lines(x_values_sim, y_values*length(x_values), lwd = 2,col = "green")
    lines(ppp[[which(multiple_p_shown=="smooth")]]$dens,lwd = 2, col = "red")
    legend('topright',
           legend=c("t","smooth"),
           lty=1,
           cex =0.5,
           trace = F,####line type
           col=c("green","red"),
           bty='n')


  }
  axis(1,pos=0)   ###the coordinate at which the axis line is to be drawn: if not NA this overrides the value of line.

  ####comment out if you do not eant axis
  #  if(side=='smaller') {axis(2,pos=0,las=1)}   #### the style of axis labels. (0=parallel, 1=all horizontal, 2=all perpendicular to axis, 3=all vertical)
#  else {axis(2,
   #          pos=h$breaks[1],las=1)}

  lines(rep(actual,2),     ###x1,x2 for the line
        c(0,h2))          ##y1 ,y2 forthe line


if(show_p==T){
  for(i in 1:length(ppp)){
  if(!is.nan(ppp[[i]]$p)&is.numeric(ppp[[i]]$p)){
    text(actual,    ###x position of the text
         h2-i*0.1*h2,        ##y position of the text
         pos=pos,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         ##respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
         labels=paste(multiple_p_shown[i],
                      ' p=',
                      signif(ppp[[i]]$p,round_number),
                      sep=''))

  }
    if(!is.nan(ppp[[i]]$p)&!is.numeric(ppp[[i]]$p)){
    text(actual,    ###x position of the text
         h2-i*0.1*h2,        ##y position of the text
         pos=pos,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         ##respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
         labels=paste(multiple_p_shown[i]," p",ppp[[i]]$p,sep='')) ####what is the text
      }
  }
}
  text(actual,    ###x position of the text
       max(h$density)*0.003,        ##y position of the text
       pos=3,
       ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
       labels=paste0(signif(actual,round_number))
  )

  if(permutation_visual=="mean"){
    text(median(distribution),    ###x position of the text
         max(h$density)*0.003,        ##y position of the text
         pos=3,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         labels=paste0("mean=",
                       signif(mean(distribution),round_number)
         )
    )
  }

  if(permutation_visual=="median"){
    text(median(distribution),    ###x position of the text
         max(h$density)*0.003,        ##y position of the text
         pos=3,
         ##a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4,
         labels=paste0("median=",
                       signif(median(distribution),round_number)
         )
    )
  }

  if(permutation_visual=="none"|missing(permutation_visual)){

  }
  ###################################################################


}


}






