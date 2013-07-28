\name{traffic}
\alias{traffic}
\docType{data}
\title{Traffic dataset
}
\description{
This dataset contains the number of drivers involved in road accidents with 
fatalities in 2005, classified according to the age and sex of the drivers and 
the vehicles involved.

The categories of  vehicles driven are the following: Bicycles, mopeds, 
motorcycles, public service cars up to 9 seats, private cars, vans, trucks 
weighing less than 3500 kg, agricultural tractors, trucks and articulated 
vehicles weighing more than 3500 kg, buses and other unspecified vehicles. 

The age categories envisaged in the study are the following: 18 to 20, 
21 to 24, 25 to 34, 35 to 44, 45 to 54, 55 to 64, 65 to 74, over 74 years 
and unknown age (UA).
}
\usage{data(traffic)}
\format{
  A data frame (11x18). 
  Rows represents the kinds of vehicle driven. 
  The first 9 columns corresponds to the ages of male drivers.
  Columns  10 to 18 corresponds to the ages of female drivers.
}
\source{
Spanish Traffic Authority \cr
http://www.dgt.es/portal/es/seguridad_vial/estadistica/accidentes_30dias/datos_desagregados.do \cr
(Grupo 4 - 2005, on sheet 4.2.C).
}
\references{
Zarraga, A. & Goitisolo, B. (2009).
   Simultaneous analysis and multiple factor analysis for contingency tables: 
   Two methods for the joint study of contingency tables.
   \emph{Computational Statistics and Data Analysis},
   \bold{53}, 3171--3182.
}
\examples{
data(traffic)
traffic
}
\keyword{ datasets }
