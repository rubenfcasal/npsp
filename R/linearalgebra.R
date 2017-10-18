.DPOSV_R <- function(a, b){
  a <- as.matrix(a)
  n <- nrow(a) 
  # square matrix
  if (!is.numeric(a) || (ncol(a) != n))
    stop("'a' is not a square numeric matrix")   
  nb <- ncol(b)  
  if ( !identical(n, nrow(b)) )
    stop("arguments 'a' and 'b' must have the same 'nrow'")
  # subroutine DPOSV_R(N, NRHS, A, B, INFO)     
  # real(8) A(N, N), B(N, NRHS)
  ret <-.Fortran("DPOSV_R", n = as.integer(n), 
                 nb = as.integer(nb), 
                 a = as.double(a), 
                 b = as.double(b), 
                 info = integer(1))
  if (ret$info > 0)
    stop("The leading minor of order", ret$info ,"of 'a' is not positive definite")
  if (ret$info < 0)
    stop("Error", ret$info ,"calling LAPACK routine DPOSV")
  # dim(ret$a) <- c(n, n)  
  # dim(ret$b) <- c(n, nb)  
  return(ret)
}


.DNRM2_R <- function(x, z){
  x <- as.matrix(x)
  nx <- nrow(x)     # number of data
  nd <- ncol(x)  
  z <- as.matrix(z)
  nz <- nrow(z)
  if ( !identical(nd, ncol(z)) )
    stop("arguments 'x' and 'z' must have the same 'ncol'")
  # subroutine DNRM2_R(nd, x, nx, z, nz, dist)
  ret <-.Fortran("DNRM2_R", nd = as.integer(nd), 
                 x = as.double(t(x)), 
                 nx = as.integer(nx), 
                 z = as.double(t(z)), 
                 nz = as.integer(nz), 
                 dist = double(nz*nx))                                   
  dim(ret$dist) <- c(nz, nx)  
  return(ret$dist)
}
                                  


