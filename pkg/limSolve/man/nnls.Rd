\name{nnls}
\alias{nnls}
\title{Nonnegative Least Squares}
\description{Solves the following inverse problem:
 \deqn{\min(||Ax-b||^2)} subject to        
 \deqn{x>=0}
 \cr uses subroutine nnls (FORTRAN) from Linpack}
\usage{nnls(A, B, tol=sqrt(.Machine$double.eps), verbose=TRUE)}
\arguments{
  \item{A }{numeric matrix containing the coefficients of the equality constraints \eqn{Ax~=B}; if the columns of A have a names attribute, they will be used to label the output}
  \item{B }{numeric vector containing the right-hand side of the equality constraints}
  \item{tol }{tolerance (for singular value decomposition, equality constraints)}
  \item{verbose }{logical to print nnls error messages}
}
\value{
  a list containing:
  \item{X }{vector containing the solution of the nonnegative least squares problem.}
  \item{residualNorm }{scalar, the sum of absolute values of residuals of violated inequalities (i.e. x<0); should be zero or very small if the problem is feasible}
  \item{solutionNorm }{scalar, the value of the quadratic function at the solution, i.e. the value of \eqn{\min(||Ax-b||^2)}}
  \item{IsError }{logical, TRUE if an error occurred}
  \item{type }{the string "nnls", such that how the solution was obtained can be traced}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{A <- matrix(nr=2,nc=3,data=c(3,2,2,4,2,1))
B <- c(-4,3)
nnls(A,B)}
\references{Lawson C.L.and Hanson R.J. 1974. Solving Least Squares Problems, Prentice-Hall
\cr Lawson C.L.and Hanson R.J. 1995. Solving Least Squares Problems.
 SIAM classics in applied mathematics, Philadelphia.        (reprint of book)}
\seealso{\code{\link{ldei}}, which includes equalities}
\keyword{algebra}
\keyword{array}
\keyword{optimize}
