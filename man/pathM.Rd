\name{pathM}
\alias{pathM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
network construction for gene list
}
\description{
Based on the expression level, this function constructs a predicted pathway map
}
\usage{
pathM(data, FC = 2, FDR = 0.05, group = "med", scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the dataset for gene expression. Consisting of a matrix. Rownames are gene names, colnames are sample ID.
}
  \item{FC}{
The fold change for difference between target gene and other genes.
}
  \item{FDR}{
The False discovery rate for results. This also called adjusted P-value by the method of Benjami and
Hochberg.
}
  \item{group}{
The group for the target gene expresion profile. default is median and group="med"; it also can set "quar" means  gene expression value less than a quartile and more than three quartile.
}
  \item{scale}{
a logical value indicating whether the expression profile needs to normalized by log2(X+1).default scale=TRUE.
}
}

\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)
dat=data[1:50,]
pathM(dat)
}

\keyword{analysis}
