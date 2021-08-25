#' Sampling of a vector
#'
#' Random sampling of a vector into n groups
#' Bug: Sampling is faulty when length(vect) < n
#' @param vect A vector to be sampled into groups
#' @param n Number of groups
#' @param sampLen A vector with custom number of samples per group. Is calculated if missing (best choice).
#' @return a list with n groups containing sub sampled `vect`
#' @export
vectSamp=function(vect,n=4,sampLen) { # Pick 'n' random samples within vector 'vect'
  # sampLen is a vector of number of observations within sample
  # If sampLen is not specified it is automatically calculated
  j=length(vect)
  if (missing(sampLen)) { # Calculation of sampLen vector
    sampLen=numeric(n)
    oddAdd=numeric(n)
    if (j%%n!=0) oddAdd[1:(j%%n)]=1    ##This is to generate remainder number of 1 ##4%%3=1,4%%5=4,4%%2=0
    sampLen=floor(j/n)+oddAdd          ##calculate number of samples per group
  }
  if (length(sampLen)!=n)
    stop("Mismatch number of samples")
  if (sum(sampLen)>j)
    stop("Exceeds maximum samples")
  if (sum(sampLen)<j)
    cat('\n Warning: Undersampling \n')
  vectSamp=sample(vect) # Randomise vector  ##takes a sample of the specified size from the elements of x using either without replacement.
  samp=list()
  vectInd=1 # Index in randomised vector
  for (i in 1:n) { # Divide randomised vector into samples
    samp[[i]]=sort(vectSamp[vectInd:(vectInd+sampLen[i]-1)])  # Make sample
    vectInd=vectInd+sampLen[i]  # Count up index
  }
  return(samp)
}
