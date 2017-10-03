#' Predict pls
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#' @param object 
#' @param newdata 
#' @param ... 
#'
#' @return pls prediction
#' @export
predict.plsMUVR <-function(object, newdata,  ...){
  #-- validation des arguments --#
  if (missing(newdata)) stop("No new data available.")
  
  x = object$X
  y = object$Y
  q = ncol(y)
  p = ncol(x)
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p) 
  }
  
  #-- initialisation des matrices --#	
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  
  means.x = attr(x, "scaled:center")
  means.y = attr(y, "scaled:center")
  sigma.x = attr(x, "scaled:scale")
  sigma.y = attr(y, "scaled:scale")
  
  newdata = as.matrix(newdata)
  ##- coeff de regression 
  B.hat = array(0, dim = c(p, q, ncomp))
  ##- prediction
  y.hat = y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
  ##- variates
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  
  variates.x = object$variates$X
  betay = list()
  
  #-- prediction --#
  for(h in 1:ncomp){
    
    dd= coefficients(lm(y~variates.x[,1:h,drop=FALSE])) #regression of y on variates.global.x => =loadings.global.y at a scale factor
    if (q==1) betay[[h]]=(dd[-1]) else betay[[h]]=(dd[-1,])
    
    W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
    B = W %*% drop(betay[[h]])
    
    y.temp=scale(newdata,center=means.x,scale=sigma.x) %*% as.matrix(B) #so far: gives a prediction of y centered and scaled
    y.temp=scale(y.temp,center=FALSE,scale=1/sigma.y) #so far: gives a prediction of y centered, with the right scaling
    y.temp=scale(y.temp,center=-means.y,scale=FALSE) #so far: gives a prediction of y with the right centering and scaling
    
    y.hat[, , h] = y.temp # we add the variance and the mean of y used in object to predict
    t.pred[, h] = scale(newdata, center = means.x, scale = sigma.x) %*% W[, h]
    B.hat[, , h] = B
  }  #end h
  
  #-- valeurs sortantes --#
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(y.hat) = rownames(newdata)
  colnames(y.hat) = colnames(y)
  
  return(invisible(list(predict = y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}


#' Predict plsda
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#' @param object 
#' @param newdata 
#' @param method 
#' @param ... 
#'
#' @return plsda predictions
#' @export
predict.plsdaMUVR=function(object, newdata, method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), ...)  {
  #-- validation des arguments --#
  if (missing(newdata)) stop("No new data available.")
  
  x = object$X
  y = object$Y 
  yprim = object$ind.mat   
  q = ncol(yprim)          
  p = ncol(x)
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ", p, 
           " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p) 
  }
  
  #-- initialisation des matrices --#	
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  
  means.x = attr(x, "scaled:center")
  means.y = attr(y, "scaled:center")
  sigma.x = attr(x, "scaled:scale")
  sigma.y = attr(y, "scaled:scale")
  
  newdata = as.matrix(newdata)
  ##- coeff de regression 
  B.hat = array(0, dim = c(p, q, ncomp))
  ##- prediction
  y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  ##- variates
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  variates.x = object$variates$X
  betay = list()
  
  #-- prediction --#
  for(h in 1:ncomp){
    dd= coefficients(lm(y~variates.x[,1:h,drop=FALSE])) #regression of y on variates.global.x => =loadings.global.y at a scale factor
    if(q==1) betay[[h]]=(dd[-1]) else betay[[h]]=(dd[-1,])
    
    W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
    B = W %*% drop(betay[[h]])
    
    y.temp=scale(newdata,center=means.x,scale=sigma.x) %*% as.matrix(B) #so far: gives a prediction of y centered and scaled
    y.temp=scale(y.temp,center=FALSE,scale=1/sigma.y) #so far: gives a prediction of y centered, with the right scaling
    y.temp=scale(y.temp,center=-means.y,scale=FALSE) #so far: gives a prediction of y with the right centering and scaling
    
    y.hat[, , h] = y.temp # we add the variance and the mean of y used in object to predict
    t.pred[, h] = scale(newdata, center = means.x, scale = sigma.x) %*% W[, h]
    B.hat[, , h] = B
  }  #end h
  
  G = matrix(0, nrow = q, ncol = ncomp)
  cls = list()
  
  for (i in 1:q) {
    if(ncomp > 1) {
      G[i, ] = apply(object$variates$X[yprim[, i] == 1, , drop = FALSE], 2, mean)
    } else {
      G[i, ] = mean(object$variates$X[yprim[, i] == 1, ])
    }
  }	
  
  # ----    max distance -----------------
  
  if (any(method == "all") || any(method == "max.dist")) {
    
    function.pred = function(x){
      nr = nrow(x)
      tmp = vector("numeric", nr)
      for(j in 1:nr){
        tmp[j] = (which(x[j, ] == max(x[j, ]))[1])
      }
      return(tmp)
    }
    cls$max.dist = matrix(apply(y.hat, 3, function.pred), ncol = ncomp)
    colnames(cls$max.dist) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
  }
  
  # ----    centroids distance -----------------
  
  if (any(method == "all") || any(method == "centroids.dist")) {
    
    cl = matrix(nrow = nrow(newdata), ncol = ncomp)
    
    centroids.fun = function(x, G, h) {
      q = nrow(G)
      x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
      if (h > 1) {
        d = apply((x - G[, 1:h])^2, 1, sum)
      }
      else {
        d = (x - G[, 1])^2
      }
      cl.id = which.min(d)
    }
    
    for (h in 1:ncomp) {
      cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, centroids.fun, G = G, h = h)
      cl[, h] = cl.id		
    }
    colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
    cls$centroids.dist = cl
  }	
  
  # ----    mahalanobis distance -----------------
  
  if (any(method == "all") || any(method == "mahalanobis.dist")) {
    
    cl = matrix(nrow = nrow(newdata), ncol = ncomp)
    
    Sr.fun = function(x, G, yprim, h) {
      q = nrow(G)
      Xe = yprim %*% G[, 1:h]
      Xr = object$variates$X[, 1:h] - Xe
      Sr = t(Xr) %*% Xr / nrow(y)
      Sr.inv = solve(Sr)
      x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
      if (h > 1) {
        mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
        d = apply(mat^2, 1, sum)
      }
      else {
        d = drop(Sr.inv) * (x - G[, 1])^2
      }
      cl.id = which.min(d)
    }
    
    for (h in 1:ncomp) {
      cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, Sr.fun, G = G, yprim = yprim, h = h)
      cl[, h] = cl.id		
    }
    colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
    cls$mahalanobis.dist = cl
  }
  
  #-- valeurs sortantes --#
  if (any(method == "all")) method = "all"
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(y.hat) = rownames(newdata)
  colnames(y.hat) = colnames(y)
  colnames(G) = paste("dim", c(1:ncomp), sep = " ")
  
  return(invisible(list(predict = y.hat, variates = t.pred, B.hat = B.hat, 
                        centroids = G, method = method, class = cls)))
}
