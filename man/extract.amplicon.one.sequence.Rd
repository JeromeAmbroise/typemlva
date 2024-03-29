% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.amplicon.R
\name{extract.amplicon.one.sequence}
\alias{extract.amplicon.one.sequence}
\title{extract.amplicon.one.sequence enables to extract from a target sequence (a DNAstringset of lenght 1) the amplicon defined by the primer pair}
\usage{
extract.amplicon.one.sequence(target, Fprobe, Rprobe, max.mismatch)
}
\arguments{
\item{target}{the reference sequence that you want to screen}

\item{Fprobe}{the sequence of the forward probe}

\item{Rprobe}{the sequence of the reverse probe}

\item{max.mismatch}{the maximum number of mismatches allowed in the probes}
}
\value{
the sequence of the amplicon (including both primers)
}
\description{
You provide both a genome and the probes
}
