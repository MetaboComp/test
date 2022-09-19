#' Stratified sampling for classification and unique individuals
#'
#' @param Y Classes (Std Y variable)
#' @param ID Individual identifier (unique)
#' @param folds Number of folds
#'
#' @return A master list of segments
#' @export
#'
#' @examples
#' Y <- rep(LETTERS[1:2],10)
#' ID <- 1:20
#' folds <- 5
#' uniqDASamp(Y, ID, folds)
uniqDASamp <- function(Y, ID, folds) {
  Ynames <- sort(unique(Y))  # Find groups
  groups <- length(Ynames)
  groupList <- list()
  for (grp in 1:groups) {
    groupID <- ID[Y==Ynames[grp]]  # Find indices per group
    groupList[[grp]] <- vectSamp(groupID,n=folds)  # Draw random samples within group
  }
  masterList <- groupList[[1]] # Add 1st groups to 'Master' sample of all groups
  for (grp in 2:groups) {  # Add subsequent groups
    masterList <- masterList[order(sapply(masterList,length))]
    for (segment in 1:folds) {
      masterList[[segment]] <- sort(c(masterList[[segment]],groupList[[grp]][[segment]]))
    }
  }
  return(masterList)
}

