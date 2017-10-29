\name{delta.estimate}
\alias{delta.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the treatment effect, the difference in survival at time t
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the treatment effect in the survival setting i.e. the difference in survival at time t between the treatment group and the control group. The inverse probability of censoring weighted estimate of survival within each treatment group is used; there is an option to use the Kaplan-Meier estimate instead. This function is generally not expected to be used directly by the user, it is called by the recover.B function.
}
\usage{
delta.estimate(xone, xzero, deltaone, deltazero, t, weight = NULL, KM = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xone}{
%%     ~~Describe \code{xone} here~~
numeric vector, the observed event times in the treatment group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{xzero}{
%%     ~~Describe \code{xzero} here~~
numeric vector, the observed event times in the control group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{deltaone}{
%%     ~~Describe \code{deltaone} here~~
numeric vector, the event indicators for the treatment group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.

}
  \item{deltazero}{
%%     ~~Describe \code{deltazero} here~~
numeric vector, the event indicators for the control group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{t}{
%%     ~~Describe \code{t} here~~
the time of interest.
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
a \eqn{n_1+n_0} by \eqn{x} matrix of weights where \eqn{n_1 =} sample size in the treatment group and \eqn{n_0 =} sample size in the control group, default is null; generally not supplied by user, only used by other functions.
}
  \item{KM}{
%%     ~~Describe \code{KM} here~~
true or false, indicating whether the Kaplan-Meier estimate of survival should be used instead of the inverse probability of censoring weighted estimate
}
}
\value{
the difference in survival at time t (treatment group minus control group)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
%%  ~~who you are~~
Layla Parast
}

\examples{
data(dataA)
delta.estimate(xone = dataA$x1, xzero = dataA$x0, deltaone = dataA$delta1, deltazero = 
dataA$delta0, t=1)

delta.estimate(xone = dataA$x1, xzero = dataA$x0, deltaone = dataA$delta1, deltazero = 
dataA$delta0, t=0.5)

}

\keyword{univar}
\keyword{nonparametric}
\keyword{survival}