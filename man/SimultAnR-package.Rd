\name{SimultAnR_1.1-package}
\alias{SimultAnR_1.1-package}
\alias{SimultAnR_1.1}
\docType{package}
\title{Correspondence and Simultaneous Analysis
} 
\description{This package performs classical correspondence analysis (CA) and 
     simultaneous analysis (SA). Simultaneous analysis is a factorial 
     methodology developed for the joint treatment of a set of several 
     contingency tables. In SA tables having the same rows, are concatenated row-wise.
     In this version of the package a multiple option have been included for the simultaneous analysis
     of tables having the same columns, concatenated column-wise.
     This way, a MSA allows to perform the analysis of an indicator matrix where the rows represent individuals. 
     In this package, functions for computation, summaries 
     and graphical visualization in two dimensions are provided, including 
     options to display partial and supplementary points.
}
\details{
\tabular{ll}{
Package: \tab SimultAnR_1.1\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-28\cr
License: \tab GPL (>=2) \cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
The main function of the package is SimAn() which performs simultaneous 
analysis.
The input dataset is the concatenation of the frequency/contingency tables to 
be jointly analyzed. The package can also be used to perform classical
correspondence analysis of frequency/contingency tables by means of the 
function CorrAn().
}
\author{Amaya Zarraga, Beatriz Goitisolo

Maintainer: Beatriz Goitisolo <beatriz.goitisolo@ehu.es>
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

Zarraga, A. & Goitisolo, B. (2011).
   Simultaneous Analysis in S-PLUS: The SimultAn Package.
   \emph{Journal of Statistical Software},
   \bold{70 (11)}, 1--22.
   
}
\keyword{ package }
\seealso{
\code{\link{CorrAn}}, 
\code{\link{summary.CorrAn}}, 
\code{\link{plot.CorrAn}}, 
\code{\link{SimAn}}, 
\code{\link{summary.SimAn}}, 
\code{\link{plot.SimAn}},    

\code{\link{shoplifting}}, 
\code{\link{traffic}}.
%%~ Optional links to other man pages, e.g. ~
%%~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~
}
\examples{
######################################################################
### Package and data                                               ###
######################################################################

data(shoplifting)

######################################################################
### Correspondence Analysis                                        ###
######################################################################
dataCA <- shoplifting[, 1:9]

### CA without supplementary elements
CorrAn.out <- CorrAn(data=dataCA)

### CA with supplementary rows/columns
CorrAn.out <- CorrAn(data=dataCA, sr=13)
CorrAn.out <- CorrAn(data=dataCA, sc=9)
CorrAn.out <- CorrAn(data=dataCA, sr=13, sc=9)

### CA with supplementary columns and without output for rows
CorrAn.out <- CorrAn(data = dataCA, sc=9, oac=0, oar = 0)

### Summary
summary(CorrAn.out)

### Graphs on screen
\donttest{plot(CorrAn.out, s1=1, s2=2, screen=TRUE)}

### Graphs on a pdf file
pdf('CAGr.pdf', paper="a4r", width=12, height=9)
plot(CorrAn.out, s1=1, s2=2, screen=FALSE)
dev.off()

###other
args(CorrAn)
names(CorrAn.out)
CorrAn.out
CorrAn.out$totalin
CorrAn.out$Gs
             
             

######################################################################
### Simultaneous Analysis                                          ###
######################################################################

dataSA <- shoplifting

### SA without supplementary elements
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:9,10:18), weight= 2,
   nameg=c("M", "F"))

### SA with supplementary rows/columns
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:8,10:17), weight= 2,
   nameg=c("M", "F"), sr= 13)
SimAn.out <- SimAn(data=dataSA, G=2, acg=list(1:8,10:17), weight= 2,
   nameg=c("M", "F"), sr= 13, sc=c(9,18))
                            
## Multiple SA with supplementary rows and without output for columns)
SimAn.out <- SimAn(data = t(dataSA), G = 2, arg = list(1:9, 10:18), 
   weight = 2, nameg=c("M", "F"), sr = 11, oac=0, multiple = 1)


### Summary
summary(SimAn.out)
                                                        
### Graphs on screen
\donttest{plot(SimAn.out)}


### Graphs on a pdf file
pdf('SAGr.pdf', paper="a4r", width=12, height=9)
plot(SimAn.out, s1=1, s2=2, screen=FALSE)
dev.off()

###other
args(SimAn)
names(SimAn.out)
SimAn.out
SimAn.out$I
SimAn.out$resi

}
\keyword{ multivariate }
\keyword{ correspondence analysis }
\keyword{ simultaneous analysis }
