% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{signifGxE}
\alias{signifGxE}
\title{Compute Z score and p-value for parameter of focus at particular
moderator level}
\usage{
signifGxE(result, focus, level)
}
\arguments{
\item{result}{output from \link{loadResults}}

\item{focus}{parameter name on which to calculate a Z score and p-value}

\item{level}{numeric level of the moderator}
}
\value{
result with new ME, ME.SE, Z, and P columns
}
\description{
Compute Z score and p-value for parameter of focus at particular
moderator level
}
\examples{
# TODO
}
\seealso{
Other reporting: 
\code{\link{isSuspicious}()},
\code{\link{loadResults}()},
\code{\link{plot.gwsemResult}()},
\code{\link{signif}()}
}
\concept{reporting}
