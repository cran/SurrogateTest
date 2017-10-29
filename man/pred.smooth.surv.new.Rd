\name{pred.smooth.surv.new}
\alias{pred.smooth.surv.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the conditional probability of survival to time t

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; calculates the estimated conditional probability of survival to time t given survival to the landmark time and given surrogate marker information. 
}
\usage{
pred.smooth.surv.new(Axzero.f, Adeltazero.f, Aszero.f, Bsnew.f, 
Bsnew2.f = NULL, myt, bw = NULL, weight.pred, extrapolate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Axzero.f}{
%%     ~~Describe \code{Axzero} here~~
observed event times in the control group in  Study A
}
  \item{Adeltazero.f}{
%%     ~~Describe \code{Adeltazero} here~~
event/censoring indicators in the control group in Study A
}
  \item{Aszero.f}{
%%     ~~Describe \code{Aszero} here~~
surrogate marker values in the control group in Study A, NA for individuals not observable at the time the surrogate marker was measured
}
  \item{Bsnew.f}{
%%     ~~Describe \code{Bsnew} here~~
surrogate marker values in one group (either treatment or control) in Study B
}
  \item{Bsnew2.f}{
%%     ~~Describe \code{Bsnew2} here~~
surrogate marker values in one group (either treatment or control) in Study B, optional
}
  \item{myt}{
%%     ~~Describe \code{myt} here~~
time of interest
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
bandwidth, optional
}
  \item{weight.pred}{
%%     ~~Describe \code{weight.pred} here~~
weight used for perturbation resampling
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Details are included in the documentation for early.delta.test
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
conditional probability of survival past t 
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{survival}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{smooth}
\keyword{robust}
\keyword{internal}
