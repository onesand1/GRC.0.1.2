\name{readEXP}
\alias{readEXP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read the data for expression profile
}
\description{
read the data for expression profile.
}
\usage{
readEXP(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The path and file name for inputing file. File type was .csv.
}
}

\value{
\item{gene }{The gene for the first column.}
\item{path}{The input file path.}
\item{data}{The expression matrix data.}
\item{sample}{The sample name for first row.}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ analysis }% use one of  RShowDoc("KEYWORDS")
