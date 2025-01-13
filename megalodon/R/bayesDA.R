bayesDA <- function(tidy_otu, groups)
{
  
  
  otu2 <- otu %>% select(Sample,Reads,NCBI,Grupo,Dieta,Tratamiento) %>%
    group_by(Sample) %>%
    mutate(TotalReads = sum(Reads), NSpecies = sum(Reads>0)) %>%
    ungroup() %>%
    group_by(NCBI) %>%
    mutate(Nzeros = sum(Reads ==0)) %>%
    ungroup() %>%
    mutate(NSamples = n_distinct(Sample)) %>%
    mutate(Freq = Reads/TotalReads) %>%
    mutate(ZeroFreq = Nzeros/NSamples) %>%
    mutate(Freq = Freq + 1e-8)
  
  
  library(doParallel)
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  
  otu3 <- otu2 %>% filter(ZeroFreq < 0.8) %>%
    select(NCBI,Sample,Freq,Grupo,Tratamiento,Dieta) %>%
    group_by(NCBI) %>%
    nest()
  
  totalRes =list()
  
  otu_tmp = otu3
  i=0
  while(nrow(otu_tmp)>1)
  {
    i = i+1
    if(nrow(otu_tmp)> 80)
    {
      tmp <- otu_tmp[1:80,]
      otu_tmp <-otu_tmp[-(1:80),]
    }else{
      tmp <- otu_tmp[1:nrow(otu_tmp),]
      otu_tmp <- data.frame()
    }
    
    
    res <- foreach (i = 1:nrow(tmp)) %dopar% {
      rstanarm::stan_betareg(Freq ~ Dieta+Grupo+Tratamiento, data = as.data.frame(tmp$data[i]))
    }
    names(res) = as.character(tmp$NCBI)
    totalRes <- append(totalRes,res)
    print(glue::glue("Ejecutado bloque {i}"))
  }
  
  stopCluster(cl)
  
  prueba <- tibble(NCBI = names(totalRes), model = totalRes) %>% inner_join(otu3)
  
  prueba2 <- prueba %>%
    mutate(test = map_dfr(model,~estimate_contrasts(model = .,contrast = c("Dieta"))))
  
  
  
}