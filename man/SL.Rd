\name{SL}
\alias{SL}

\title{
Diferential analysis of mRNA expression
}
\description{
Diferential analysis of mRNA expression
}
\usage{
SL(data,group = "group",exp_type = "rpkm",plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{data}{
a data frame,including 3 or more columns ('sample','gene name','group') for one gene or many genes.
}
  \item{group}{
a character of group column name for sample groups
}

  \item{exp_type}{
the expression data formula mode. including 'RPKM','FPKM','TPM'.
}
  \item{plot}{
a logical value for ploting results.default TRUE, print the figure to "DE_plot" directory.}
}

\value{
\item{gene }{Gene names or IDs.}
\item{P }{P value for differential analysis based on expression.}
\item{FC }{Fold change for differential analysis based on expression.}
\item{FC_type }{The Fold change type.}
\item{method }{The differential analysis method.}

}


\examples{
data=system.file("extdata","data.Rda", package="GRC")
load(data)
data=t(data)
data1=data[,1:2]
group=as.numeric(substr(rownames(data),14,15))
group[which(group>9)]="Normal"
group[which(group<10)]="Cancer"
sample=rownames(data)
dat=cbind(sample,data1)
dat=cbind(dat,group)

colnames(dat)=c('sample',colnames(data1),'group')
SL(data=dat,plot = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")
