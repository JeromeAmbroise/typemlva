#' extract.amplicon.one.sequence enables to extract from a target sequence (a DNAstringset of lenght 1) the amplicon defined by the primer pair
#'
#' You provide both a genome and the probes
#'
#' @param target the reference sequence that you want to screen
#' @param Fprobe the sequence of the forward probe
#' @param Rprobe the sequence of the reverse probe
#' @param max.mismatch the maximum number of mismatches allowed in the probes
#'
#' @return the sequence of the amplicon (including both primers)
#' @import Biostrings
#'
#' @export
extract.amplicon.one.sequence <- function(target,Fprobe,Rprobe,max.mismatch)
{
  library(Biostrings)
  start <- NULL
  end <- NULL

  match <- vmatchPattern(pattern =Fprobe,subject = target,max.mismatch = max.mismatch)
  start <- c(start,start(unlist(match)))
  end <- c(end,end(unlist(match)))

  match <- vmatchPattern(pattern =as.character(reverseComplement(DNAString(Fprobe))),subject = target,max.mismatch = max.mismatch)
  start <- c(start,start(unlist(match)))
  end <- c(end,end(unlist(match)))

  match <- vmatchPattern(pattern =Rprobe,subject = target,max.mismatch = max.mismatch)
  start <- c(start,start(unlist(match)))
  end <- c(end,end(unlist(match)))

  match <- vmatchPattern(pattern =as.character(reverseComplement(DNAString(Rprobe))),subject = target,max.mismatch = max.mismatch)
  start <- c(start,start(unlist(match)))
  end <- c(end,end(unlist(match)))

  start.end <- c(start,end)

  if(length(start.end)==4)
  {
    minimum <- min(c(start.end))
    maximum <- max(c(start.end))

    subsequence <- subseq(target,start=minimum,end = maximum)

    return(subsequence)
  }

}




#' extract.amplicon.multiple.sequence enables to extract from a target sequence (a DNAstringset of lenght n) the amplicon(s) defined by the primer pair
#'
#' You provide both a genome and the probes
#'
#' @param target the reference sequence that you want to screen
#' @param Fprobe the sequence of the forward probe
#' @param Rprobe the sequence of the reverse probe
#' @param max.mismatch the maximum number of mismatches allowed in the probes
#'
#' @return the sequence of the amplicon (including both primers)
#' @import Biostrings
#'
#' @export
extract.amplicon.multiple.sequence <- function(target,Fprobe,Rprobe,max.mismatch)
{
  library(Biostrings)
  subsequence <- DNAStringSet()
  for(i in 1:length(target))
  {subsequence <- c(subsequence,extract.amplicon.one.sequence(target=target[i],Fprobe='TGACTACTGAAACAGTTTTTG',Rprobe='ATGATTGTACCGAGTAAAAGA',max.mismatch=1))}
  return(subsequence)
}
