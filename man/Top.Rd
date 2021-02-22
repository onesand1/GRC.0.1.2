\name{Top}
\alias{Top}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Achieve  top n genes for target
}
\description{
When the number of results produced by "regulator" function is too large, this function can extract the pre-sorted results from above results.
}
\usage{
Top(object, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The results from regulator function.
}
  \item{n}{
The number for display genes order by FC.
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

\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")

