#' Title
#'
#' @param file_list
#'
#' @return
#' @export
#'
#' @examples
import_kraken <- function(file_list, threads)
{
  require(doParallel)

  res <- data.frame() %>% as_tibble

  cl <- makeCluster(threads)
  registerDoParallel(cl)
  res <- foreach( f = file_list, .combine = rbind) %dopar% {
  {
    table <- read.delim(f, sep = "\t", header = F, stringsAsFactors = F)
    colnames(table) <- c("Prop","Reads","Assigned","TaxRank","NCBI","Name")
    table$Sample = basename(f)
    table$depth <- nchar(gsub("\\S.*","",table$Name))/2
    table$Name2 <- gsub("^ *","",table$Name)
    table$depth = table$depth+1
    table$Name3 <- paste0(table$TaxRank,"__",table$Name2)

    last = c()
    for(i in 1:nrow(table))
    {
      last[table$depth[i]] = table$Name2[i]
      table$Fullname[i] = paste(last[1:table$depth[i]],sep = ";",collapse = ";")
    }
    last = c()
    for(i in 1:nrow(table))
    {
      last[table$depth[i]] = table$Name3[i]
      table$Fullname2[i] = paste(last[1:table$depth[i]],sep = ";",collapse = ";")
    }


    last = c()
    for(i in 1:nrow(table))
    {
      last[table$depth[i]] = table$NCBI[i]
      table$FullNCBI[i] = paste(last[1:table$depth[i]],sep = ";",collapse = ";")
    }
    res <- rbind(res,table)
  }
  }
  res <- res %>% distinct()
  stopCluster(cl)
  return(res)
}
