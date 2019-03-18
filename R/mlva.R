type.mlva <- function(genome,start,end,motif)
{
  nhit <- numeric()

  for(i in 1:30)
  {
    sequence.to.search <- paste0(start,paste(rep(motif,i),collapse = ''),end)
    mycount <- vcountPattern(pattern = sequence.to.search,subject = genome,max.mismatch = 0)
    nhit[i] <- sum(mycount)
  }

  result <- which(nhit==1)
  return(result)
}
