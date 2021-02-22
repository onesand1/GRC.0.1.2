\name{cutoff}
\alias{cutoff}
\title{
cutoff grouping for survival analysis
}
\description{
Desgin three method to identify the cutoff for survival analysis, including median,  quartile and Youden index.

}
\usage{
cutoff(data, type = "M", ttype = "d", cut = 1,plot=TRUE,file="Rplot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data needs three columns: expression value, survival times and status. The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
}
  \item{type}{
Type refers to the segmentation threshold of survival analysis, which is divided into three types: median(M), quartile (Q) and Youden index (Y).default type="M".
}
  \item{ttype}{
This parameter refers to the units of time, days(d), months(m), and years(y).default time='d'.
}
  \item{cut}{
This parameter is a threshold used when selecting type= "Y". marker values to use as a cut-off for calculation of sensitivity and specificity.
}
  \item{plot}{
A logical parameter. This parameter is defined to show whether the survival curve is visualized.default plot=TRUE.
}
  \item{file}{
The survival curve plot file name. The file type is .pdf.
}
}

\value{
\item{fit}{ a survfit object.}
 \item{data}{ Added a group column for raw data.}
%% ...
}

\examples{
##
library(survival)
library(GRC)
library(survminer)
library(survivalROC)

###
data(lung)
n=dim(lung)[1]
data=cbind(lung$age,lung$time)
data=cbind(data,lung$status)
cutoff(data,type = "Y")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")
