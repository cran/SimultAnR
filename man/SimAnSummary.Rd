\name{SimAnSummary}
\alias{SimAnSummary}
\title{Summary of Simultaneous Analysis
}
\description{This function summarizes the results of \code{SimAn}.
}
\usage{
SimAnSummary(input)
}
\arguments{
  \item{input }{The output of the simultaneous analysis }
}
\details{
The function \code{SimAnSummary} gives the detailed numerical results of the 
\code{SimAn} function. 

In its first stage simultaneous analysis performs a simple correspondence 
analysis of each table, so the summary contains the separate correspondence 
analysis of each table as provided by the \code{CorrAn} function. 

The joint analysis of all the tables is performed in the second stage of the 
simultaneous analysis and total inertia, the eigenvalues, percentages of 
explained inertia and cumulated percentages of explained inertia for all 
dimensions are listed. 
The output also contains for the overall rows and for the columns of the 
tables, the masses, chi-squared distances and, by default restricted to the 
first two dimensions, projections of points on each dimension or principal 
coordinates, contributions of the points to the dimensions and squared 
correlations.
For partial rows and for supplementary elements the same results are listed 
except for contributions of the points to the dimensions.

The output of \code{SimAnSummary} also contains the relations between overall 
rows and partial rows, the relations between the factors of the CA of the 
different tables, the relations between the factors of the SA and the factors 
of the separate CA of the different tables, the projections of the tables and 
the contributions of each table 
to the principal axes.
}
\value{
Results of separate correspondence analysis of each table:
  \item{   Total inertia }{
     Total inertia, as a measure of the total variance of the data table      }
  \item{   Eigenvalues and percentages of inertia }{
     Eigenvalues or principal inertias and percentages of explained inertia   }
  \item{   Output for rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension, contributions and squared correlations         }
  \item{   Output for columns }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension, contributions and squared correlations         }
  \item{   Output for supplementary rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
  \item{   Output for supplementary columns  }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
     
Results of simultaneous analysis of the set of tables:
  \item{   Total inertia }{
     Total inertia, as a measure of the total variance of the data table      }
  \item{   Eigenvalues and percentages of inertia }{
     Eigenvalues or principal inertias }
  \item{   Output for rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension, contributions and squared correlations         }
  \item{   Output for columns }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension, contributions and squared correlations         }
  \item{   Output for partial rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
  \item{   Projections of tables}{
     Projections of each table on each dimension                              }
  \item{   Contributions of tables to SA}{
     Contributions of each table to the dimensions                            }
  \item{   Relation between overall and partial rows}{
     Relation between overall and partial rows                                }
  \item{   Relation between factors of separate CA}{
     Relation between factors of separate correspondence analysis             }
  \item{   Relation between factors of CA and SA}{
     Relation between factors of correspondence analysis and simultaneous 
     analysis                                                                 }
  \item{   Output for supplementary rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
  \item{   Output for supplementary partial rows }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
  \item{   Output for supplementary columns }{
     Masses, chi-squared distances of points to their average, projections of 
     points on each dimension and squared correlations                        }
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
\code{\link{SimAnGraph}}.
}
\examples{
data(shoplifting)
dataSA <- shoplifting

### SA without supplementary elements
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:9,10:18), weight= 2,
   nameg=c("M", "F"))

### Summary
SimAnSummary(SimAn.out)
}
\keyword{ multivariate }
\keyword{ simultaneous analysis}

