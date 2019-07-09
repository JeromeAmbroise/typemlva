
#' typemlva enables to compute the mlva profile of an isolate based on the genome
#'
#' You provide both a genome and the motifs. The function computes the number of repetition for each motif
#'
#' @param genome the reference genome that you want to screen
#' @param motifs a list with each motif
#' @param kmer the kmer size which was used during genome assembly
#'
#' @return a characther with the MLVA profile
#' @import Biostrings
#'
#' @export
typemlva <- function(genome,motifs,kmer)
{
  library(Biostrings)
  genome <- genome[width(genome)>2000]
  m1 <- type1motif(genome,motif=motifs[1],kmer)
  m2 <- type1motif(genome,motif=motifs[2],kmer)
  m3 <- type1motif(genome,motif=motifs[3],kmer)
  m4 <- type1motif(genome,motif=motifs[4],kmer)
  m5 <- type1motif(genome,motif=motifs[5],kmer)
  profile <- paste0('(',paste(m1,m2,m3,m4,m5,sep=';'),')')
  return(profile)
}

#' type1motif enables to compute the number of motif for a single VNTR locus
#'
#' You provide both a genome and the motif. The function computes the number of repetition for this motif
#'
#' @param genome the reference genome that you want to screen
#' @param motifs the sequence of the motif
#' @param kmer the kmer size which was used during genome assembly
#'
#' @return the number of repetitions
#' @import Biostrings
#'
#' @export
type1motif <- function(genome,motif,kmer)
{

  nhit.forward <- numeric()
  nhit.reverse <- numeric()
  index <- 100
  i <- 1

  while(index >0)
  {
    sequence.to.search <- paste0(paste(rep(motif,i), collapse = ""))
    mycount.forward <- vcountPattern(pattern = sequence.to.search,subject = genome, max.mismatch = 0)
    nhit.forward[i] <- sum(mycount.forward)
    mycount.reverse <- vcountPattern(pattern = as.character(reverseComplement(DNAString(sequence.to.search))),subject = genome, max.mismatch = 0)
    nhit.reverse[i] <- sum(mycount.reverse)
    index <- nhit.forward[i] + nhit.reverse[i]
    i <- i + 1
  }

  nhit <- nhit.forward + nhit.reverse
  result <- max(which(nhit >= 1))
  if(result>=floor(kmer/nchar(motif)))
  { result <- floor(kmer/nchar(motif))
  result <- paste(expression('>='),as.character(result))}

  return(result)
}



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
  {subsequence <- c(subsequence,extract.amplicon(target=target[i],Fprobe='TGACTACTGAAACAGTTTTTG',Rprobe='ATGATTGTACCGAGTAAAAGA',max.mismatch=1))}
  return(subsequence)
}
