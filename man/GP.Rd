\name{GP}
\alias{GP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 the Corelated Analysis between gene and pathway
}
\description{
  The corelated analysis between genesets (pathways and GOs) and gene  based on expression profile.
}
\usage{
GP(expression, geneset,log=TRUE, gsva,plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{
a matrix for expression profile
}
  \item{geneset}{
  a string for gmt file path or a list for one or more gene lists.

}
  \item{log}{
a logical value indicating whether the expression profile needs to normalized by log2(X+1).default log=TRUE.
}
  \item{gsva}{
the option for GSVA analysis refer to the package GSVA.
}
  \item{plot}{
A logical parameter. This parameter is defined to show whether the correlation is visualized.The figure data will save to "cor_plot" directory.default plot=TRUE.
}
}

\value{
 \item{cor_res}{ The analysis results of spearman correlation between pathway/Go Score and gene expression}
  \item{geneSet.name}{ The Geneset label}
  \item{geneSet.gene}{ The gene list of Geneset}
   \item{geneSet.score}{ The Score for Geneset by GSEA algorithm}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)

gene="100134869"
data=as.matrix(data)
res=regulator(gene,data)
gg=list(TEST=res$target[c(1:12)])
GP(data,gg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}
