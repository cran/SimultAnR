\name{SimAn}
\alias{SimAn}
\title{Simultaneous Analysis
}
\description{
Simultaneous analysis is a factorial method developed for the joint treatment 
of a set of several data tables, especially frequency tables whose row margins 
are different, for example when the tables are from different samples or 
different time points, without modifying the internal structure of each table. 
In the data tables rows must refer to the same entities, but columns may be 
different.
With the multiple option tables having the same columns are concatenated 
column-wise.
This way, a MSA allows to perform the analysis of an indicator matrix 
where the rows represent individuals. 
}
\usage{
SimAn(data, G, acg, weight = 2, nameg = NA, sr = NA, sc = NA, 
nd = 2, dp = 2, oar = 1, oac = 1, multiple = 0, arg)
}
\arguments{
  \item{data    }{Data set                                           }                                   
  \item{G       }{Number of tables to be jointly analyzed            }                                   
  \item{acg     }{List of number of the active columns for each table
                  (if multiple = 0)}
  \item{weight  }{Weighting on each table                            }
  \item{nameg   }{Prefix for identifying partial rows and tables     }
  \item{sr      }{Indices of supplementary rows                      }
  \item{sc      }{Indices of supplementary columns                   }
  \item{nd      }{Number of dimensions in results                    }
  \item{dp      }{Number of digits in results                        }
  \item{oar     }{Output for active rows (1 = yes, 0 = no)           }
  \item{oac     }{Output for active columns (1 = yes, 0 = no)        }
  \item{multiple}{Multiple Simultaneous Analysis (1 = yes, 0 = no)   }
  \item{arg     }{List of number of the active rows for each table   
                 (if multiple = 1)}
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
  \item{resj    }{Results of active columns                                   }
  \item{resig   }{Results of partial rows
                  (if multiple = 0) }
  \item{resjg   }{Results of partial columns 
                  (if multiple = 1)                }
  \item{Fsg     }{Projections of each table                                   }
  \item{ctrg    }{Contribution of each table to the axes                      }
  \item{riig    }{Relation between the overall rows and the partial rows  
                  (if multiple = 0)                                           }
  \item{rjjg    }{Relation between the overall rows and the partial columns
                 (if multiple = 1)                                            }
  \item{RCACA   }{Relation between separate CA axes                           }
  \item{RCASA   }{Relation between CA axes and SA axes                        }
  \item{Fs      }{Projections of active rows                                  }
  \item{Gs      }{Projections of active columns                               }
  \item{Fsig    }{Projections of partial rows                                
                  (if multiple = 0)                                           }
  \item{Gsjg    }{Projections of partial columns                                 
                  (if multiple = 1)                                           }
  \item{allFs   }{Projections of rows and partial rows (if multiple = 0) 
                   in an array format                                         } 
  \item{allGs   }{Projections of columns and partial columns (if multiple = 1)
                   in an array format                                         }
  \item{I       }{Number of active rows                                       
                  (if multiple = 0)                                           }
  \item{J       }{Number of active columns                                       
                  (if multiple = 1)                                           }
  \item{maxJg   }{Maximum number of columns for a table                       
                  (if multiple = 0)                                           }
  \item{maxIg   }{Maximum number of rows for a table                            
                  (if multiple = 1)                                           }
  \item{G       }{Number of tables                                            }
  \item{namei   }{Names of active rows                                        
                  (if multiple = 0)                                           }
  \item{namej   }{Names of active columns                                        
                  (if multiple = 1)                                           }
  \item{nameg   }{Prefix for identifying partial points, tables, etc          }
  \item{resisr  }{Results of supplementary rows                               }
  \item{resjsc  }{Results of supplementary columns                            }   
  \item{resigsr }{Results of partial supplementary rows                       
                  (if multiple = 0)                                           }
  \item{resjgsc }{Results of partial supplementary columns                       
                  (if multiple = 1)                                           }
  \item{Fssr    }{Projections of supplementary rows                           }
  \item{Gssc    }{Projections of supplementary columns                        }
  \item{Fsigsr  }{Projections of partial supplementary rows                   
                  (if multiple = 0)                                           }
  \item{Gsjgsc  }{Projections of partial supplementary columns                   
                  (if multiple = 1)                                           }
  \item{allFssr }{Projections of supplementary rows and partial supplementary 
                  rows (if multiple = 0)  in an array format                  }                 
  \item{allGssc }{Projections of supplementary columns and partial 
                  supplementary columns (if multiple = 1) in an array format  }    
  \item{Isr     }{Number of supplementary rows                                
                  (if multiple = 0)                                           }
  \item{Jsc     }{Number of supplementary columns                             
                  (if multiple = 1)                                           }
  \item{nameisr }{Names of supplementary rows                                 
                  (if multiple = 0)                                           }
  \item{namejsc }{Names of supplementary columns                              
                  (if multiple = 1)                                           }
  \item{CAres   }{Results of CA of each table to be used in Summary and Graph 
                  functions }
  \item{multiple}{Value of option multiple}
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

Zarraga, A. & Goitisolo, B. (2011).
   Simultaneous Analysis in S-PLUS: The SimultAn Package.
   \emph{Journal of Statistical Software},
   \bold{70 (11)}, 1--22.
}

\author{Amaya Zarraga, Beatriz Goitisolo
}


\seealso{
\code{\link{summary.SimAn}}, 
\code{\link{plot.SimAn}}. 
}
\examples{
data(shoplifting)
dataSA <- shoplifting


### SA without supplementary elements
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:9,10:18), weight= 2, 
   nameg=c("M", "F"))

### Multiple SA without output for columns
SimAn.out <- SimAn(data=t(dataSA), G=2, weight= 2, 
   nameg=c("M", "F"), oac=0, multiple=1, arg=list(1:9,10:18))

### Summary
summary(SimAn.out)

### Graphs on screen
\donttest{plot(SimAn.out)}

### Graphs on a pdf file (without columns)
pdf('SAGr.pdf', paper="a4r", width=12, height=9)
plot(SimAn.out, s1=1, s2=2, screen=FALSE, oac=0)
dev.off()
}
\keyword{ multivariate }
\keyword{ simultaneous analysis }

