#' Plot MV object
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#' @param factCols An optional vector with colors for the factor levels (in the same order as the levels)
#' @param sampLabels Sample labels (optional; implemented for classification)
#' @param ylim Optional for imposing y-limits for regression and classification analysis
#'
#' @return A plot of results from multivariate predictions
#' @export
plotMV <- function(MVObj,
                   model = 'mid',
                   factCols,
                   sampLabels,
                   ylim = NULL) {
  
  # Basic sanity
  if (!any(class(MVObj) == 'MVObject')) stop('Wrong object class')
  
  # Model number
  modNum <- ifelse(model == 'min',
                   1,
                   ifelse(model == 'mid',
                          2,
                          3))
  # Extract actual data
  Y <- MVObj$inData$Y
  nSamp <- length(Y)
  if(missing(sampLabels)) sampLabels <- Y
  
  # Sanity check sample labels
  if(length(sampLabels)!=nSamp) stop('Length of sampLabels not equal to number of samples in Y.')

  if (any(class(MVObj) == 'Regression')) {

    ###########################
    # REGRESSION PLOT
    ###########################

    # Y-predicted overall
    YP <- MVObj$yPred[,modNum]
    # Y-predicted per repetition
    YPR <- MVObj$yPredPerRep[[modNum]]
    # Y-limits
    if (is.null(ylim)) ylim <- range(YPR)
    # Plot Y-predicted per repetition in grey
    matplot(Y,
            YPR,
            pch = 20,
            xlab = 'Original Y',
            ylab = 'Predicted Y',
            col = 'grey',
            bty = 'l',
            cex = 0.5, 
            ylim = ylim)
    # Add in overall Y-predictions in black
    points(Y,
           YP,
           pch = 20)
    # Add simple regression line
    reg <- lm(YP ~ Y)
    abline(reg)
    # Add legend
    legend('topleft',
           legend = c(paste('Model R2 =', signif(MVObj$fitMetric$R2[modNum], 3)),
                      paste('Model Q2 =', signif(MVObj$fitMetric$Q2[modNum], 3))),
           bty='n')
  } else if (any(class(MVObj) == 'Classification')) {
    
    ################################
    # CLASSIFICATION SWIMLANE PLOT
    ################################
    
    # Y-predicted overall
    YP <- MVObj$yPred[[modNum]]
    # Y-predicted per repetition
    YPR <- MVObj$yPredPerRep[[modNum]]
    # Y-limits
    if (is.null(ylim)) ylim <- range(YPR)
    # Unique levels in Y
    classes <- 1:length(levels(Y))
    # Colors per level
    if(missing(factCols)) factCols=classes+1
    if(length(factCols)!=length(classes)) stop('Length of factCols not equal to number of levels in Y.')
    # Sort out "jitter"/nudge between levels in the swimlane plot
    classNudge <- 0.2 * ((classes - mean(classes)) / (mean(classes) - 1))
    # Allocate plot surface
    plot(1:nSamp,
         Y,
         type = 'n',
         ylim = ylim,
         xlab = '',
         ylab = 'Class prediction probability',
         xaxt = 'n')
    # Custom axis
    axis(1,
         at = 1:length(Y),
         labels = sampLabels,
         las = 3)
    # Plot each Y level separately
    for(cl in classes) {
      # Y-pred per rep
      matpoints((1:nSamp) + classNudge[cl],
                YPR[,cl,],
                pch = 20,
                col = factCols[cl],
                cex = 0.5)
      # Y-pred overall
      points((1:nSamp) + classNudge[cl],
             YP[,cl],
             pch = 20,
             col = factCols[cl])
    }
    # Add swimlane lines
    for (li in 1:(nSamp + 1)) {
      abline(v = li - .5,
             lty = 3,
             col = 'grey')
    }
    # Identify erroneous classifications
    yClass <- MVObj$yClass[,modNum]
    whichWrong <- which(yClass!=Y)
    wrongClass <- as.numeric(Y[whichWrong])
    # Mark them out in the plot
    for (w in 1:length(wrongClass)) {
      points(whichWrong[w] + classNudge[wrongClass[w]],
             YP[whichWrong[w], wrongClass[w]],
             cex = 2)
    }
    # Add legend
    xpdOld <- par()$xpd
    par(xpd = TRUE)
    legend(x = 0,
           y = ylim[2] + diff(ylim) / 5,
           horiz = TRUE,
           legend = c(levels(Y),'misclassified'),
           pch = c(rep(16,length(classes)), 1),
           col = c(factCols,1),
           cex = 0.8,
           pt.cex = c(rep(0.5,length(classes)),2),
           bty = 'n')
    par(xpd = xpdOld)
  } else if (any(class(MVObj) == 'Multilevel')){
    
    ###########################
    # MULTILEVEL PLOT
    ###########################
    
    # Y-predicted overall
    YP <- MVObj$yPred[,modNum]
    # Y-predicted per repetition
    YPR <- MVObj$yPredPerRep[[modNum]]
    # Plot Y-predicted per repetition in grey
    matplot(YPR,
            1:nSamp,
            pch = 20,
            col = 'grey',
            cex = 0.5,
            ylim = c(nSamp, 1),
            ylab = 'Sample number',
            xlab = 'Predicted Y')
    # Plot Y-predicted overall in black
    points(YP,
           1:nSamp,
           pch = 20,
           col = 'black')
    # Draw support lines
    abline(h = nSamp / 2 + 0.5, lty = 2)
    abline(v = 0, lty = 2)
  }
}
