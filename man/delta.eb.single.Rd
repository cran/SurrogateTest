\name{delta.eb.single}
\alias{delta.eb.single}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the early treatment effect estimate in Study B
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the early treatment effect estimate in Study B, generally not called directly by the user
}
\usage{
delta.eb.single(Axzero, Adeltazero, Aszero, Bxzero, Bdeltazero, 
Bszero, Bxone, Bdeltaone, Bsone, t, landmark, weightA = NULL, 
weightB = NULL, weight.both = NULL, extrapolate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Axzero}{
%%     ~~Describe \code{Axzero} here~~
observed event times in the control group in  Study A
}
  \item{Adeltazero}{
%%     ~~Describe \code{Adeltazero} here~~
event/censoring indicators in the control group in Study A
}
  \item{Aszero}{
%%     ~~Describe \code{Aszero} here~~
surrogate marker values in the control group in Study A, NA for individuals not observable at the time the surrogate marker was measured
}
  \item{Bxzero}{
%%     ~~Describe \code{Bxzero} here~~
observed event times in the control group in  Study B
}
  \item{Bdeltazero}{
%%     ~~Describe \code{Bdeltazero} here~~
event/censoring indicators in the control group in Study B
}
  \item{Bszero}{
%%     ~~Describe \code{Bszero} here~~
surrogate marker values in the control group in Study B, NA for individuals not observable at the time the surrogate marker was measured
}
  \item{Bxone}{
%%     ~~Describe \code{Bxone} here~~
observed event times in the treatment group in  Study B
}
  \item{Bdeltaone}{
%%     ~~Describe \code{Bdeltaone} here~~
event/censoring indicators in the treatment group in  Study B
}
  \item{Bsone}{
%%     ~~Describe \code{Bsone} here~~
surrogate marker values in the treatment group in Study B, NA for individuals not observable at the time the surrogate marker was measured

}
  \item{t}{
%%     ~~Describe \code{t} here~~
time of interest
}
  \item{landmark}{
%%     ~~Describe \code{landmark} here~~
landmark time of interest, t0
}
  \item{weightA}{
%%     ~~Describe \code{weightA} here~~
Study A weights used for perturbation resampling
}
  \item{weightB}{
%%     ~~Describe \code{weightB} here~~
Study B weights used for perturbation resampling

}
  \item{weight.both}{
%%     ~~Describe \code{weight.both} here~~
Study A and Study B weights used for perturbation resampling
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
early treatment effect estimate
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(dataA)
data(dataB)

delta.eb.single(Axzero = dataA$x0, Adeltazero = dataA$delta0, 
Aszero = dataA$s0, Bxzero = dataB$x0, Bdeltazero = dataB$delta0, 
Bszero = dataB$s0, Bxone = dataB$x1, Bdeltaone = dataB$delta1, 
Bsone = dataB$s1, t=1, landmark=0.5,  extrapolate = TRUE)

}


\keyword{survival}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{internal}
