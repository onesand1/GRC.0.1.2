\name{regulator}
\alias{regulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlations analysis of single gene
}
\description{
Correlations analysis of single gene
}
\usage{
regulator(target, data, plot = TRUE, type = "1", scale = TRUE,
group = "med", FC = 2, FDR = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the gene that will be check in expression profile.
}
  \item{data}{
the dataset for gene expression. Consisting of a matrix. Rownames are gene names, colnames are sample ID.
}
  \item{plot}{
a logical value for ploting results.default TRUE, print the figure to window.
}
  \item{type}{
The ragulation type for target. if type=1,  thie model will caculate the downregulated genes for target gene; if type=2, it will caculate the  upregulated genes for target gene.
}
  \item{scale}{
a logical value indicating whether the expression profile needs to normalized by log2(X+1).default scale=TRUE.
}
  \item{group}{
The group for the target gene expresion profile. default is median and group="med"; it also can set "quar" means  gene expression value less than a quartile and more than three quartile.
}
  \item{FC}{
The fold change for difference between target gene and other genes.
}
  \item{FDR}{
The False discovery rate for results.if type=2,the FDR equal to P-value.
}
}


\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)
gene="100134869"
data=as.matrix(data)
res=regulator(gene,data)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")
