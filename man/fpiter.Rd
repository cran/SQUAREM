\name{fpiter}
\alias{fpiter}
\title{Fixed-Point Iteration Scheme}
\description{A function to implement the fixed-point iteration algorithm.  This includes monotone, contraction mappings including EM and MM algorithms}
\usage{
  fpiter(par, fixptfn, objfn=NULL, control=list( ), ...)
}
\arguments{
    \item{par}{A vector of parameters denoting the initial guess for the
       fixed-point iteration.}

    \item{fixptfn}{A vector function, \eqn{F}{F} that denotes the fixed-point 
       mapping. This function is the most essential input in the package.  It 
       should accept a parameter vector as input and should return a parameter
       vector of same length. This function defines the fixed-point iteration: 
       \eqn{x_{k+1} = F(x_k)}{x[k+1] = F(x[k]}.  
       In the case of EM algorithm, \eqn{F}{F} defines a single E and M step.} 

    \item{objfn}{This is a scalar function, $L$, that denotes a ''merit'' 
       function which attains its local minimum at the fixed-point of $F$.  
       This function should accept a parameter vector as input and should 
       return a scalar value. In the EM algorithm, the merit function \eqn{L}{L}
       is the log-likelihood. In some problems, a natural merit function may 
       not exist, in which case the algorithm works with only \code{fixptfn}.
       The merit function function \code{objfn} does not have to be specified,
       even when a natural merit function is available, especially when its 
       computation is expensive.} 

    \item{control}{ A list of control parameters to pass on to the algorithm.  Full names of control list elements must be specified, otherwise, user-specifications are ignored.  See *Details* below.}

    \item{...}{Arguments passed to \code{fixptfn} and \code{objfn}.}
}

\value{A list with the following components:
\item{par}{Parameter,\eqn{x*}{x*} that are the fixed-point of \eqn{F}{F} such 
   that \eqn{x* = F(x*)}{x* = F(x*)}, if convergence is successful.}

\item{value.objfn}{The value of the objective function \eqn{L}{L} at termination.}

\item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}

\item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.  
  \code{0} indicates successful convergence, 
  whereas \code{1} denotes failure to converge.}
}

\details{
\code{control} is list of control parameters for the algorithm.  

\describe{
\code{control = list(tol = 1.e-07, maxiter = 1500, trace = FALSE)}

    \code{tol}{  A small, positive scalar that determines when iterations 
       should be terminated.  Iteration is terminated when 
       \eqn{||x_k - F(x_k)|| \leq tol}{abs(x[k] - F(x[k]) <= tol}.  
       Default is \code{1.e-07}.}

    \code{maxiter}{  An integer denoting the maximum limit on the number of 
       evaluations of  \code{fixptfn}, \eqn{F}{F}.  Default is \code{1500}.}

    \code{trace}{  A logical variable denoting whether some of the intermediate 
       results of iterations should be displayed to the user.  
       Default is \code{FALSE}.} 
  }
}

\seealso{
  \code{\link{squarem}}
   }
   
\examples{

##############################################################################
# Example 1:  EM algorithm for Poisson mixture estimation 
poissmix.em <- function(p,y) {
# The fixed point mapping giving a single E and M step of the EM algorithm
# 
pnew <- rep(NA,3)
i <- 0:(length(y)-1)
zi <- p[1]*exp(-p[2])*p[2]^i / (p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i)
pnew[1] <- sum(y*zi)/sum(y)
pnew[2] <- sum(y*i*zi)/sum(y*zi)
pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
p <- pnew
return(pnew)
}

poissmix.loglik <- function(p,y) {
# Objective function whose local minimum is a fixed point \
# negative log-likelihood of binary poisson mixture
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return ( -sum(loglik) )
}

# Real data from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))
y <- poissmix.dat$freq
tol <- 1.e-08

# Use a preset seed so the example is reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=54321))

p0 <- c(runif(1),runif(2,0,4))  # random starting value

# Basic EM algorithm
pf1 <- fpiter(p=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))


##############################################################################
# Example 2:  Accelerating the convergence of power method iteration for 
# finding the dominant eigenvector of a matrix 

power.method <- function(x, A) {

# Defines one iteration of the power method
# x = starting guess for dominant eigenvector
# A = a square matrix

ax <- as.numeric(A \%*\% x)
f <- ax / sqrt(as.numeric(crossprod(ax)))
f
}

# Finding the dominant eigenvector of the Bodewig matrix
# This is a famous matrix for which power method has trouble converging
# See, for example, Sidi, Ford, and Smith (SIAM Review, 1988) 
#
# Note: there are two eigenvalues that are equally dominant, 
#  but have opposite signs.
# Sometimes the power method finds the eigenvector corresponding to the 
# large positive eigenvalue, but other times it finds the eigenvector
# corresponding to the large negative eigenvalue
b <- c(2, 1, 3, 4, 1,  -3,   1,   5,  3,   1,   6,  -2,  4,   5,  -2,  -1)
bodewig.mat <- matrix(b,4,4)
eigen(bodewig.mat)

p0 <- rnorm(4)

# Standard power method iteration
ans1 <- fpiter(p0, fixptfn=power.method, A=bodewig.mat)
# re-scaling the eigenvector so that it has unit length
ans1$par <- ans1$par / sqrt(sum(ans1$par^2))  
ans1

}

\references{ 
R Varadhan and C Roland (2008), Simple and globally convergent numerical 
   schemes for accelerating the convergence of any EM algorithm, 
   \emph{Scandinavian Journal of Statistics}, 35:335-353.

C Roland, R Varadhan, and CE Frangakis (2007), Squared polynomial extrapolation 
   methods with cycling: an application to the positron emission tomography 
   problem, \emph{Numerical Algorithms}, 44:159-172.
}

\keyword{optimization, EM algorithm}

