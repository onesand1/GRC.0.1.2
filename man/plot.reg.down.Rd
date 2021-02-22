\name{plot.reg.down}
\alias{plot.reg.down}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visual function for down-regulated results
}
\description{
Visual function for down-regulated results.
}
\usage{
plot.reg.down(object, r = 2, type = "line")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
the regulator function,GL or GP result object.
}
  \item{r}{
radius for every gene backgroud.
}
  \item{type}{
the figure type. line represent the regulated genes in line; circle represent the regulated genes in circle.
}
}

\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)
gene="100134869"
data=as.matrix(data)
res=regulator(gene,data)
plot.reg.down(res,type="circle")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}
