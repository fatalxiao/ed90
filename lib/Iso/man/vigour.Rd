\name{vigour}
\alias{vigour}
\docType{data}
\title{
  vigour
}
\description{
Growth vigour of stands of spruce trees in New Brunswick, Canada.
}
\usage{data("vigour")}
\format{
  A data frame with 23 observations (rows).  The first column
  is the year of observation (1965 to 1987 inclusive).  The other
  five columns are observations on the vigour of growth of the
  given stand in each of the years.
}
\details{
  The stands each had different initial tree densities.  It was
  expected that vigour would initially increase (as the trees
  increased in size) and then level off and start to decrease as the
  growing trees encroached upon each others' space and competed more
  strongly for resources such as moisture, nutrients, and light.
  It was further expected that the position of the mode of the
  vigour observations would depend upon the initial densities.
}
  
\source{
  These data were collected and generously made available by Kirk
  Schmidt who was at the time of collecting the data a graduate
  student in the Department of Forest Engineering at the University
  of New Brunswick, Fredericton, New Brunswick, Canada.  The data
  were collected as part of his research for his Master's degree
  (supervised by Professor Ted Needham) at the University of New
  Brunswick.  See Schmidt (1993).
}
\references{
K. D. Schmidt  (1993). \emph{Development of a precommercial
thinning guide for black spruce}.  Thesis (M.Sc.F.), University of
New Brunswick, Faculty of Forestry.
}
\examples{
 matplot(vigour[,1],vigour[,2:6],
         main="Growth vigour of stands of New Brunswick spruce",
         xlab="year",ylab="vigour",type="b")
}
\keyword{datasets}
