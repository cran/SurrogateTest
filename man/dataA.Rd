\name{dataA}
\alias{dataA}
\docType{data}
\title{
Hypothetical Study A data
}
\description{
Hypothetical Study A data to be used in examples; t=1 and the landmark time = 0.50.

}
\usage{data(dataA)}
\format{
 A list with 6 elements representing 1000 observations from a control group and 1000 observations from a treatment group:
  \describe{
    \item{\code{s1}}{Surrogate marker measurement for treated observations; this marker is measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, this value is NA.}
    \item{\code{x1}}{The observed event or censoring time for treated observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta1}}{The indicator identifying whether the treated observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
   \item{\code{s0}}{Surrogate marker measurement for control observations; this marker is measured at time = 0.5. For observations that experience the primary outcome or are censored before 0.5, this value is NA.}
    \item{\code{x0}}{The observed event or censoring time for control observations; X  = min(T, C) where T is the time of the primary outcome and C is the censoring time.}
    \item{\code{delta0}}{The indicator identifying whether the control observation was observed to have the event or was censored; D =1*(T<C) where T is the time of the primary outcome and C is the censoring time.}
}}
\details{
Note that if the observation is censored or experienced the primary outcome before the landmark time of 0.50, the surrogate marker measurement is not observed and coded NA. 
}
\examples{
data(dataA)
names(dataA)
}
\keyword{datasets}
