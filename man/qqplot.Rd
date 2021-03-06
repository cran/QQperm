% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QQ.R
\name{qqplot}
\alias{qqplot}
\title{QQ plot of observed P-values vs expected P-values, using the empirical (permutation-based) expected p-value distribution. This empirical-based expected p-value distribution no longer depends on an assumption that the Fisher's Exact two-tailed p-values are uniformly distributed under the null. For a given matrix, the permutation-based expected distribution is plotted relative to the observed order statistic to get the permutation-based QQ plot.}
\usage{
qqplot(P.perm, P.observed, adjust.xy = TRUE, ...)
}
\arguments{
\item{P.perm}{Expected P-values from NULL distribution, which is generated through permutation in our example.}

\item{P.observed}{Observed P-values from true case/control assignment.}

\item{adjust.xy}{An option to have the x-axis and y-axis adjusted based on their own range in the plot.}

\item{...}{Additional arguments passed to the plot function.}
}
\value{
None.
}
\description{
QQ plot of observed P-values vs expected P-values, using the empirical (permutation-based) expected p-value distribution. This empirical-based expected p-value distribution no longer depends on an assumption that the Fisher's Exact two-tailed p-values are uniformly distributed under the null. For a given matrix, the permutation-based expected distribution is plotted relative to the observed order statistic to get the permutation-based QQ plot.
}
\examples{
library(QQperm)
data("example.Ps")

#print output to pdf file only if not running in interactive mode
if (!interactive()) {
 pdf("qqplot.pdf")
}

qqplot(example.Ps$perm, example.Ps$observed)

if (!interactive()) {
 dev.off()
}
}
\author{
Slave Petrovski and Quanli Wang
}

