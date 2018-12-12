###-----------------------------------------------------------------------------
### Projective Decomposition:
###
### The scale-invariant matrix Z equivalent to the (m x n) matrix A
### is defined by the (projective decomposition) equation:
###        A     = s Diag{r} Z Diag{c}
###     { a_ij } = { s * r_i * z_ij * c_j }
### where
###     s = rms(A),
###     r_i > 0, rms(z_i*) = 1 for all i = 1..m, and
###     c_j > 0, rms(z_*j) = 1 for all j = 1..n
### Here rms(x) = sqrt(mean(x*x)) denotes the root-mean-square of all elements
### in a vector or matrix.
###
### Projective Decomposition defines an equivalence relation among matrices,
### in the sense that the scale-invariant matrix equivalent to a matrix is
### unique (if it exists), and any two matrices X and Y with the same
### scale-invariant matrix and the same pattern of signs are equivalent
### matrices, i.e. there are strictly positive vectors a of length m and b of
### length n so that Diag{a} X Diag{b} = Y.
###-----------------------------------------------------------------------------
eps <- 1e-7

# Coefficient of variation; the data MUST all be nonnegative
cv  <- function(x) { sqrt(var(x))/max(mean(x),eps) }

# Root-mean-square: Euclidean distance normalized for the number of dimensions
rms <- function(x) { sqrt(mean(x*x)) }

###
# E[A] = projective.scale(projective.decomposition(A)) = s r t(c)
#
# This is the expected value of A under a model in which the data is
# centered at the origin; variation around the origin is modeled as a
# scaled outer product of the vector of relative row scaling
# factors (r) and the vector of relative column scaling factors (c).
###
projective.scale <- function(D)
{
  D$rms * D$row.factor %*% t(D$col.factor)
}

###
# Z = scale.invariant(projective.decomposition(A)) = A / (s r t(c))
#
# This is the expected value of A under a model in which the data is
# centered at the origin; variation around the origin is modeled as a
# scaled outer product of the vector of relative row scaling
# factors (r) and the vector of relative column scaling factors (c).
###
scale.invariant <- function(A,D)
{
  A / projective.scale(D)
}

###
# Comparing the maximum of two (non-negative) values to a tolerance
# ensures both values are within the tolerance
#     row.radii <- sqrt(pmax(rowMeans(M*M,na.rm=TRUE),0)) # rms of each row
#     col.radii <- sqrt(pmax(colMeans(M*M,na.rm=TRUE),0)) # rms of each column
#     row.rv    <- cv(row.radii) # Variability of the row radii, relative to their average
#     col.rv    <- cv(col.radii) # Variability of the column radii, relative to their average
#     asphericity <- max(row.rv,col.rv)
###
asphericity <- function(M) {
  max(cv(sqrt(pmax(rowMeans(M*M,na.rm=TRUE),0))),
      cv(sqrt(pmax(colMeans(M*M,na.rm=TRUE),0))))
}

###
# D <- projective.decomposition(A)
#
# This function computes the _model_, i.e. the parameters necessary
# to convert between A and its normalized version Z.
#
# The meanings of the returned components, in terms of the projective
# decomposition equation
#     A = s Diag{r} Z Diag{c}
# are:
#     D$rms        = s = rms(A)
#     D$row.factor = r
#     D$col.factor = c
# These total m+n+1 values, which is typically much less than the
# size of A or Z (m*n). When Z is required, rather than just the
# parameters, one can either compute
#    Z = scale.invariant(A,projective.decomposition(A)),
# or compute each z_ij = s r_i a_ij c_j as it is required.
###
projective.decomposition <- function(A,
                                     tol      = 1e-6,
                                     max.iter = 100,
                                     verbose  = FALSE,
                                     digits   = 3)
{
  m <- dim(A)[1]
  n <- dim(A)[2]
  rms.A <- rms(A)
  B <- A/rms.A
  r <- rep(1.0,m)
  c <- rep(1.0,n)
  for (iter in c(1:max.iter)) {
    W    <- B / (r[row(B)]*c[col(B)])
    asph <- asphericity(W)
    if (verbose) {
      print(paste(iter,signif(asph/tol,digits)))
    }
    if (asph < tol) { break }
    r <- r * sqrt(pmax(rowMeans(W*W,na.rm=TRUE),eps))
    W <- B / (r[row(B)]*c[col(B)])
    c <- c * sqrt(pmax(colMeans(W*W,na.rm=TRUE),eps))
    g <- sqrt((rms(c)+eps) / (rms(r)+eps))
    r <- r * g
    c <- c / g
  }
  # return value
  list('rms'        = rms.A,
       'row.factor' = r,
       'col.factor' = c,
       'iterations' = iter,
       'tolerance'  = asph
  )
}

# end of Projective.R
