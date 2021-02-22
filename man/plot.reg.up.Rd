\name{plot.reg.up}
\alias{plot.reg.up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visual function for up-regulated results
}
\description{
Visual function for up-regulated results
}
\usage{
plot.reg.up(object, r = 2, type = "line")
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
the figure type. "line" represent the regulated genes in line; "circle" represent the regulated genes in circle.
}
}

\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)
gene="100134869"
data=as.matrix(data)
res=regulator(gene,data,type=2)
res=Top(res,10)
plot.reg.up(res,type="circle")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")
