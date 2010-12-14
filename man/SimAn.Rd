\name{SimAn}
\alias{SimAn}
\title{Simultaneous Analysis
}
\description{
Simultaneous analysis is a factorial method developed for the joint treatment 
of a set of several data tables, especially frequency tables whose row margins 
are different, for example when the tables are from different samples or 
different time points, without modifying the internal structure of each table. 
In the data tables  rows must refer to the same entities, but columns may be 
different.
}
\usage{
SimAn(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA, 
nd = 2, dp = 2)
}
\arguments{
  \item{data    }{Data set                                           }                                   
  \item{G       }{Number of tables to be jointly analyzed            }                                   
  \item{acg     }{List of number of the active columns for each table}
  \item{weight  }{Weighting on each table                            }
  \item{nameg   }{Prefix for identifying partial rows and tables     }
  \item{sr      }{Indices of supplementary rows                      }
  \item{sc      }{Indices of supplementary columns                   }
  \item{nd      }{Number of dimensions in results                    }
  \item{dp      }{Number of digits in results                        }
}
\details{
The parameter \code{weight} refers to the weighting of each table included in 
simultaneous analysis in order to balance the influence of each table in the 
joint analysis, as measured by the inertia, and to prevent the joint analysis 
from being dominated by a particular table. The choice of this weighting 
depends on the aims of the analysis and on the initial structure of the 
information, and different values may be used. 
Three values are possible, \code{weight = 1} means no weighting , 
\code{weight = 2} means that the weighting is the inverse of the first 
eigenvalue (square of first singular value) of each table and is given by 
default, and \code{weight = 3} means that the weighting is the inverse of the 
total inertia of each table.

The parameter \code{nameg} allows the user to distinguish in the interpretation 
of the results as well as in the graphical representations which partial rows 
belong to each table.
By default, if this parameter is not indicated, partial rows of the first table 
will  be identified as \code{G1} followed by the name of the row, partial rows 
of the second table as \code{G2} followed by the name of the row and so on. 
The \code{nameg} argument also allows the different tables in the analysis to 
be identified.
}
\value{
  \item{totalin }{Total inertia                                               }
  \item{resin   }{Results of inertia                                          }
  \item{resi    }{Results of active rows                                      }
  \item{resig   }{Results of partial rows                                     }
  \item{resj    }{Results of active columns                                   }
  \item{Fsg     }{Projections of each table                                   }
  \item{ctrg    }{Contribution of each table to the axes                      }
  \item{riig    }{Relation between the overall rows and the partial rows      }
  \item{RCACA   }{Relation between separate CA axes                           }
  \item{RCASA   }{Relation between CA axes and SA axes                        }
  \item{Fsi     }{Projections of active rows                                  }
  \item{Fsig    }{Projections of partial rows                                 }
  \item{Gs      }{Projections of active columns                               }
  \item{allFs   }{Projections of rows and partial rows in an array format     }
  \item{allGs   }{Projections of columns in an array format                   }
  \item{I       }{Number of active rows                                       }
  \item{maxJg   }{Maximum number of columns for a table                       }
  \item{G       }{Number of tables                                            }
  \item{namei   }{Names of active rows                                        }
  \item{nameg   }{Prefix for identifying partial rows, tables, etc            }
  \item{resisr  }{Results of supplementary rows                               }
  \item{resigsr }{Results of partial supplementary rows                       }
  \item{Fsisr   }{Projections of supplementary rows                           }
  \item{Fsigsr  }{Projections of partial supplementary rows                   }
  \item{allFssr }{Projections of rows and partial supplementary rows in an 
                  array format                                                }
  \item{Isr     }{Number of supplementary rows                                }
  \item{nameisr }{Names of supplementary rows                                 }
  \item{resjsc  }{Results of supplementary columns                            }
  \item{Gssc    }{Projections of supplementary columns                        }
  \item{allGssc }{Projections of supplementary columns in an array format     }
  \item{Jsc     }{Number of supplementary columns                             }
  \item{namejsc }{Names of supplementary columns                              }
  \item{CAres   }{Results of CA of each table to be used in Summary and Graph 
                  functions }
}

\references{
Goitisolo, B. (2002).
   \emph{El Analisis Simultaneo. Propuesta y aplicacion de un nuevo metodo de 
   analisis factorial de tablas de contingencia}. Phd thesis, Basque Country 
   University Press, Bilbao.

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

\author{Amaya Zarraga, Beatriz Goitisolo
}


\seealso{
\code{\link{SimAnSummary}}, 
\code{\link{SimAnGraph}}. 
}
\examples{
data(shoplifting)
dataSA <- shoplifting

### SA without supplementary elements
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:9,10:18), weight= 2, 
   nameg=c("M", "F"))

### SA with supplementary rows/columns
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:8,10:17), weight= 2,
   nameg=c("M", "F"), sr= 13)
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:8,10:17), weight= 2,
   nameg=c("M", "F"), sr= 13, sc=c(9,18))

### Summary
SimAnSummary(SimAn.out)

### Graphs on screen
\donttest{SimAnGraph(SimAn.out)}

### Graphs on a pdf file
pdf('SAGr.pdf', paper="a4r", width=12, height=9)
SimAnGraph(SimAn.out, s1=1, s2=2, screen=FALSE)
dev.off()
}
\keyword{ multivariate }
\keyword{ simultaneous analysis }

