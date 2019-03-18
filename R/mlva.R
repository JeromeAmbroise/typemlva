typemlva <- function (genome, start, end, motif)
{
  nhit.forward <- numeric()
  nhit.reverse <- numeric()


  for (i in 1:30) {
    sequence.to.search <- paste0(start, paste(rep(motif,i), collapse = ""), end)
    mycount.forward <- vcountPattern(pattern = sequence.to.search,subject = genome, max.mismatch = 1)
    nhit.forward[i] <- sum(mycount.forward)
    mycount.reverse <- vcountPattern(pattern = as.character(reverseComplement(DNAString(sequence.to.search))),subject = genome, max.mismatch = 1)
    nhit.reverse[i] <- sum(mycount.reverse)
  }

  nhit <- nhit.forward + nhit.reverse

  result <- which(nhit == 1)
  return(result)
}
