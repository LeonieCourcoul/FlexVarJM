% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.manag.long.R
\name{data.manag.long}
\alias{data.manag.long}
\title{Management of longitudinal data}
\usage{
data.manag.long(formGroup, formFixed, formRandom, data.long1)
}
\arguments{
\item{formGroup}{A formula which indicates the group variable}

\item{formFixed}{A formula which indicates the fixed effects for the longitudinal submodel}

\item{formRandom}{A formula which indicates the random effects for the longitudinal submodel}

\item{data.long1}{A dataframe with the longitudinal data}
}
\value{
A list with the following components :
\describe{
\item{\code{data_long}}{a clean dataframe for the longitudinal data}
\item{\code{y.new.prog}}{the vector of responses variable}
\item{\code{X}}{a matrix with the fixed effects}
\item{\code{U}}{a matrix with the random effects}
\item{\code{id}}{a vector with the identification of individuals}
\item{\code{offset}}{a vector with the number of measurements for each individual}
\item{\code{I}}{an integer, the number of individuals}
}
}
\description{
Management of longitudinal data
}
