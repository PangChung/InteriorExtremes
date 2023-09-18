\name{tmvnsim}
\alias{tmvnsim}
\title{Truncated Multivariate Normal Simulation}
\description{
  This function simulates an importance sample from the truncated multivariate 
  normal distribution with mean equal to \code{mean} and covariance matrix
  \code{sigma}, lower and upper truncation points \code{lower} and \code{upper}.
}

\usage{
tmvnsim(nsamp, k, lower=rep(-Inf,k), upper=rep(Inf,k), imod=rep(FALSE, k)
		, means=rep(0, k), sigma=diag(1, k))
tmvtsim(nsamp, k, lower=rep(-Inf,k), upper=rep(Inf,k), imod=rep(FALSE, k)
		, means=rep(0, k), sigma=diag(1, k),df=nu)
}

\arguments{
 \item{nsamp}{Number of samples to generate}
 \item{k}{Dimension}
 \item{lower}{Vector of lower truncation points, 
              default is \code{rep(-Inf, length = k)}.}
 \item{upper}{Vector of upper truncation points, 
              default is \code{rep( Inf, length = k)}.}
 \item{imod}{Logical vector indicating whether truncation bounds of a 
             particular dimension should be applied on the modulus (absolute value) 
	     of that component. Defaults to FALSE for all components}
 \item{means}{Mean vector, default is 0 for all components}
 \item{sigma}{Covariance matrix, default is identity matrix}
 \item{df}{Degrees of freedom for the t-distribution, default is 4}
}

\details{

The simulation of truncated multivariate normal is done using importance sampling
(GHK simulator) that uses a Cholesky decomposition of the covariance matrix to 
convert the problem of covering the truncation region to a recursive sampling from
certain univariate truncated normals. Importance sampling weights are simultaneously 
returned. The original GHK algorithm has been slightly modified to allow truncation 
of some or all components using absolute value (modulus).
}
\value{
A list with the following two components:
\item{samp}{A matrix of nsamp rows and k columns containing the simulated samples}
\item{wts}{A vector of length nsamp containing the importance sampling weights}
}
\examples{
low <- rep(1, 5)
high <- rep(2, 5)
sig <- matrix(0.1, 5, 5)
diag(sig) <- 1
res <- tmvnsim(nsamp=1000, k=5, lower=low, upper=high, imod=rep(TRUE, 5)
		, means=rep(0, 5), sigma=sig)
head(res$samp)
head(res$wts)
cat("Probability of the region:\n")
mean(res$wts)
cat("Mean of Z-bar conditonal on truncation region:\n")
sum(rowMeans(res$samp) * res$wts)/sum(res$wts)
}
\author{Samsiddhi Bhattacharjee <sb1@nibmg.ac.in>}
\references{
Geweke, J. 1989. Bayesian inference in econometric models using Monte Carlo 
integration. Econometrica 57: 1317-1339.

Hajivassiliou, V., and D. McFadden. 1998. The method of simulated scores for the
estimation of LDV models. Econometrica 66: 863-896.

Keane, M. P. 1994. A computationally practical simulation estimator for panel data.
Econometrica 62: 95-116
}

