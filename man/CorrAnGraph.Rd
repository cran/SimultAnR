\name{CorrAnGraph}
\alias{CorrAnGraph}
\title{Graph of Correspondence Analysis 
}
\description{Graphical display of correspondence analysis results in two 
dimensions.
}
\usage{
CorrAnGraph(input, s1 = 1, s2 = 2, screen = TRUE)
}
\arguments{
  \item{input  }{The output of the correspondence analysis}
  \item{s1     }{Dimension to plot on horizontal axis}
  \item{s2     }{Dimension to plot on vertical axis}
  \item{screen }{TRUE (all the graphs on screen) or FALSE (graphs in a file)}
}
\details{
The parameters \code{s1 = 1} and \code{s2 = 2} indicate that the graph is given 
for the first two dimensions.  
A plot of e.g. the second and the third dimensions is obtained by setting 
\code{s1 = 2} and \code{s2 = 3}. 
Notice that in this case the argument \code{nd} of the \code{CorrAn} function 
must be at least 3.
If there are supplementary elements in the analysis, the \code{CorrAnGraph} 
function provides two graphical outputs, one for the active elements and the 
second one for both active and supplementary elements.
}
\author{Amaya Zarraga, Beatriz Goitisolo
}

\seealso{
\code{\link{CorrAn}}, 
\code{\link{CorrAnSummary}}.
}
\examples{
data(shoplifting)
dataCA <- shoplifting[, 1:9]

CorrAn.out <- CorrAn(data=dataCA)

### Graphs on screen
CorrAnGraph(CorrAn.out, s1=1, s2=2, screen=TRUE)

### Graph on a pdf file
pdf('CAGr.pdf',paper="a4r",width=12, height=9)
CorrAnGraph(CorrAn.out, s1=1, s2=2, screen=FALSE)
dev.off()

}
\keyword{ multivariate }
\keyword{ correspondence analysis }

