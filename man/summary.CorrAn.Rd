\name{summary.CorrAn}
\alias{summary.CorrAn}
\title{Summary of Correspondence Analysis
}
\description{This function summarizes the results of \code{CorrAn}.
}
\usage{
\method{summary}{CorrAn}(object, oar = 1, oac = 1, ...)
}
\arguments{
  \item{object}{The output of the correspondence analysis (class 'CorrAn')}
  \item{oar   }{Output for active rows (1 = yes, 0 = no)}
  \item{oac   }{Output for active columns (1 = yes, 0 = no)}
  \item{...   }{Further arguments passed to or from other methods}
}
\details{
The function \code{summary.CorrAn} gives the detailed numerical results 
of the \code{CorrAn} function corresponding to the total inertia, as a 
measure of the total variance of the data table, to the 
eigenvalues or principal inertias as well as to the  
percentages of explained inertia  and cumulated percentages of 
explained inertia for all possible dimensions. The output also 
contains, for rows and columns, the masses in \%, the 
chi-squared distances of points to their average 
and, by default restricted to the first two dimensions,  
the projections of points on each dimension or principal 
coordinates, contributions of the points to the dimensions and squared 
correlations. 
}
\value{
  \item{Total inertia}{
        Total inertia, as a measure of the total variance of the data table}
  \item{Eigenvalues and percentages of inertia }{
        Eigenvalues or principal inertias and percentages of explained inertia}
  \item{Output for rows }{
        Masses, chi-squared distances of points to their average, projections 
        of points on each dimension, contributions and squared correlations}
  \item{Output for columns }{
        Masses, chi-squared distances of points to their average, projections 
        of points on each dimension, contributions and squared correlations}
  \item{Output for supplementary rows }{
        Masses, chi-squared distances of points to their average, projections 
        of points on each dimension and squared correlations}
  \item{Output for supplementary columns }{
        Masses, chi-squared distances of points to their average, projections 
        of points on each dimension and squared correlations}
}
\references{
Greenacre, M. (2007). 
   \emph{Correspondence Analysis in Practice}. 2nd edition. Chapman and 
   Hall/CRC, London.

Lebart, L; Piron, M., Morineau, A. (2006).  
   \emph{Statistique exploratoire multidimensionnelle: visualisations et 
   inferences en fouille de donnees}. 4th edition. Dunod, Paris.
}
\author{Amaya Zarraga, Beatriz Goitisolo
}

\seealso{
\code{\link{CorrAn}}, 
\code{\link{plot.CorrAn}}.
}
\examples{
data(shoplifting)
dataCA <- shoplifting[, 1:9]

### CA without supplementary elements
CorrAn.out <- CorrAn(data=dataCA)

### Summary without output for rows
summary(CorrAn.out, oar=0)
}
\keyword{ multivariate }
\keyword{ correspondence analysis }

