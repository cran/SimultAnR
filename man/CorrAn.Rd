\name{CorrAn}
\alias{CorrAn}
\title{Correspondence Analysis
}
\description{The \code{CorrAn} function computes the correspondence analysis of 
the selected data.
}
\usage{
CorrAn(data, sr = NA, sc = NA, nd = 2, dp = 2)
}
\arguments{
  \item{data }{Data set                         }
  \item{sr   }{Indices of supplementary rows    }
  \item{sc   }{Indices of supplementary columns }
  \item{nd   }{Number of dimensions in results  }
  \item{dp   }{Number of digits in results      }
}
\details{
The options \code{sr} and \code{sc} allow supplementary rows and columns to be 
specified.
}
\value{
  \item{totalin }{Total inertia                                                }
  \item{eig     }{Eigenvalues                                                  }
  \item{resin   }{Results of inertia                                           }
  \item{resi    }{Results of active rows                                       }
  \item{resj    }{Results of active columns                                    }
  \item{resisr  }{Results of supplementary rows                                }
  \item{resjsc  }{Results of supplementary columns                             }
  \item{X       }{Matrix X                                                     }
  \item{totalk  }{Total of data table                                          }
  \item{I       }{Number of active rows                                        }
  \item{namei   }{Names of active rows                                         }
  \item{fi      }{Marginal of active rows                                      }
  \item{Fs      }{Projections of active rows                                   }
  \item{d2i     }{Chi-square distance of active rows to their average          }
  \item{J       }{Number of active columns                                     }
  \item{namej   }{Names of active columns                                      }
  \item{fj      }{Marginal of active columns                                   }
  \item{Gs      }{Projections of active columns                                }
  \item{d2j     }{Chi-square distance of active columns to their average       }
  \item{Isr     }{Number of supplementary rows                                 }
  \item{nameisr }{Names of supplementary rows                                  }
  \item{fisr    }{Marginal of supplementary rows                               }
  \item{Fssr    }{Projections of supplementary rows                            }
  \item{d2isr   }{Chi-square distance of supplementary rows to the   average   }
  \item{Xsr     }{Matrix X for supplementary rows                              }
  \item{Jsc     }{Number of supplementary columns                              }
  \item{namejsc }{Names of supplementary columns                               }
  \item{fjsc    }{Marginal of supplementary columns                            }
  \item{Gssc    }{Projections of supplementary columns                         }
  \item{d2jsc   }{Chi-square distance of supplementary columns to the   average}     
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
\code{\link{CorrAnGraph}}, 
\code{\link{CorrAnSummary}}  
}
\examples{

data(shoplifting)
dataCA <- shoplifting[1:13, 1:9]

### CA without supplementary elements
CorrAn(data=dataCA)

### CA with supplementary rows/columns
CorrAn(data=dataCA, sr=13, sc=9)

}
\keyword{ multivariate }
\keyword{ correspondence analysis }

