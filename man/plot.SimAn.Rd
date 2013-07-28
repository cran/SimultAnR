\name{plot.SimAn}
\alias{plot.SimAn}
\title{Graphs of Simultaneous Analysis
}
\description{Graphical representation of simultaneous analysis.
}
\usage{
\method{plot}{SimAn}(x, s1 = 1, s2 = 2, screen = TRUE, oar = 1, oac = 1, ...)
}
\arguments{
  \item{x      }{The output of the simultaneous analysis (class 'SimAn')}
  \item{s1     }{Dimension to plot on horizontal axis}
  \item{s2     }{Dimension to plot on vertical axis}
  \item{screen }{TRUE (all the graphs on screen) or FALSE (graphs in a file)\cr 
                 FALSE for S-plus}
  \item{oar    }{Output for active rows (1 = yes, 0 = no)}
  \item{oac    }{Output for active columns (1 = yes, 0 = no)}
  \item{...    }{Further arguments passed to or from other methods}
}
\details{
The parameters \code{s1 = 1} and \code{s2 = 2} (by default) indicate that the 
graph is given for the first two dimensions.  
A plot of e.g. the second and the third dimensions is obtained by setting 
\code{s1 = 2} and \code{s2 = 3}. 
Notice that in the latter case \code{nd}, the dimensionality of the solution, 
must be at least 3.

A large number of elements can take part in the SA of a set of tables depending 
on the number of tables and on the dimensions of the tables. 
In order to facilitate interpretation,  the \code{plot.SimAn} function gives 
several graphical outputs for active and supplementary  points (if any) 
including the ones corresponding to the separate CA of each table.

Each table is identified with a different color and the columns and partial 
rows of each table are in the same color and have the same symbol.
Two additional graphs, one for the projections of the tables and one for the 
relations between factors of CA and SA are also provided by \code{plot.SimAn}.
}
\references{
Goitisolo, B. (2002).
   \emph{El Analisis Simultaneo. Propuesta y aplicacion de un nuevo metodo de 
   analisis factorial de tablas de contingencia}. 
   Phd thesis, Basque Country University Press, Bilbao.

Zarraga, A. & Goitisolo, B. (2002).
   Methode factorielle pour l analyse simultanee de tableaux de contingence.
   \emph{Revue de Statistique Appliquee},    
   \bold{L}, 47--70

Zarraga, A. & Goitisolo, B. (2003).
   Etude de la structure inter-tableaux a travers l Analyse Simultanee,
   \emph{Revue de Statistique Appliquee},
   \bold{LI}, 39--60.

Zarraga, A. and Goitisolo, B. (2006).
   Simultaneous analysis: A joint study of several contingency tables with 
   different margins.
   In: M. Greenacre, J. Blasius (Eds.),
   \emph{Multiple Correspondence Analysis and Related Methods}, Chapman & 
   Hall/CRC, Boca Raton, Fl, 327--350.

Zarraga, A. & Goitisolo, B. (2009).
   Simultaneous analysis and multiple factor analysis for contingency tables: 
   Two methods for the joint study of contingency tables.
   \emph{Computational Statistics and Data Analysis},
   \bold{53}, 3171--3182.
}

\author{Amaya Zarraga, Beatriz Goitisolo}

\seealso{
\code{\link{SimAn}}, 
\code{\link{summary.SimAn}}.
}
\examples{
data(shoplifting)
dataSA <- shoplifting

SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:9,10:18), weight= 2,
   nameg=c("M", "F"))
                                                       
### Graphs on screen without columns
\donttest{plot(SimAn.out, oac=0)}

### Graph on a pdf file
pdf('SAGr.pdf',paper="a4r",width=12,height=9)
plot(SimAn.out, s1=1, s2=2, screen=FALSE)
dev.off()
}
\keyword{ multivariate }
\keyword{ simultaneous analysis }

