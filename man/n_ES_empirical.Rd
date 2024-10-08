% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-generation.R
\name{n_ES_empirical}
\alias{n_ES_empirical}
\title{Simulate empirical distribution of sample size and number of effect
sizes}
\usage{
n_ES_empirical(dat)
}
\arguments{
\item{dat}{a \code{data.frame} or \code{tibble} containing primary study sample
sizes and number of effect sizes per primary study}
}
\value{
A function that generates random samples from the input dataset.
}
\description{
A functional that takes in a dataset with empirical distribution of
primary study sample sizes and number of effect sizes per primary study and
returns a function that generates random samples from the dataset
}
\examples{
study_features <- n_ES_empirical(wwc_es)
study_features(m=3)
study_features(m=7)

}
