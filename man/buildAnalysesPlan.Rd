% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{buildAnalysesPlan}
\alias{buildAnalysesPlan}
\title{Build a plan for data analyses}
\usage{
buildAnalysesPlan(snpData, sliceSize)
}
\arguments{
\item{snpData}{a pathway to a file containing GWAS data. The data can be 
in a variety of forms, such as standard PLINK format (bed/bim/fam), 
PLINK2 format (pgen/pvar/psam), Oxford format (bgen/sample), or CSV 
format (csv format in much slower due to the lack of compression 
for non-binary files).}

\item{sliceSize}{number of SNPs to analyze per job}
}
\value{
Returns a data.frame with one job specification per row with the following columns:
\describe{
\item{path}{Path to the genetic data file}
\item{begin}{Starting SNP}
\item{end}{Ending SNP}
\item{count}{Number of SNPs in this job}
\item{slice}{Within data file slice index}
}
}
\description{
Long-running jobs are vulnerable to early termination from
maintanance or power outages. We recommend chopping your analyses
into smaller chunks. This also offers the advantage of running jobs
in parallel. This function builds a plan that roughly splits the
whole analysis into equal amounts of work.
}
\examples{
dir <- system.file("extdata", package = "gwsem")
buildAnalysesPlan(file.path(dir,"example.bgen"), 45)
}
