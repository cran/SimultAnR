\name{shoplifting}
\alias{shoplifting}
\docType{data}
\title{Shoplifting dataset
}
\description{
This dataset contains data from the CBS Statistics on shoplifting (Israels, 
1987). 

The data concerns persons suspected of shoplifting in 1977 and 1978 in Dutch 
stores and big textile shops, classified according to the sex and age of the 
person and the kind of stolen goods. 

The categories of stolen goods are the following: clothing, clothing accesory, 
provisions and/or tobacco, writing materials, books, records, household goods, 
sweets, toys, jewelry, perfume, hobby and/or tools and other.

The age categories envisaged in the study are the following: less than 12, 
12 to 14, 15 to 17, 18 to 20, 21 to 29, 30 to 39, 40 to 49, 50 to 64, 65 or 
over.
}
\usage{data(shoplifting)}
\format{
  A data frame (13x18). 
  Rows represent the kinds of stolen goods.
  The first 9 columns correspond to the ages of male suspected of shoplifting.
  Columns 10 to 18 correspond to the ages of female suspected of shoplifting.
}
\references{

Israels, A. (1987).
   \emph{Eigenvalue Techniques for Qualitative Data}.
   DSWO Press, Leiden.

Zarraga, A. & Goitisolo, B. (2002).
   Methode factorielle pour l analyse simultanee de tableaux de contingence.
   \emph{Revue de Statistique Appliquee},    
   \bold{L}, 47--70

Zarraga, A. & Goitisolo, B. (2003).
   Etude de la structure inter-tableaux a travers l Analyse Simultanee,
   \emph{Revue de Statistique Appliquee},
   \bold{LI}, 39--60.

}

\examples{
data(shoplifting)
shoplifting
}
\keyword{ datasets }
