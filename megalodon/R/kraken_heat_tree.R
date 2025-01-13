kraken_heat_tree <- function(deseq_results,tax_table, threshold)
{
  taxtable <- tax_table %>% select(id = NCBI, TaxRank,Name2) %>% distinct() %>% mutate(id = as.character(id))
  tax_network <- deseq_results %>%
  filter(padj < threshold) %>%
  rownames_to_column("NCBI") %>%
  inner_join(tax_table %>%
               mutate(NCBI = as.character(NCBI)) %>%
               select(NCBI,FullNCBI,Name2) %>%
               distinct())

    net <- data.frame()
  for(i in 1:nrow(tax_network))
  {
    c <- str_split(tax_network$FullNCBI[i],pattern = ";")[[1]] %>% as.vector()
    if(length(c)>1){
      for(j in 1:(length(c)-1))
      {
        net <- bind_rows(net, data.frame(from = c[j], to =c[j+1]))
      }
    }
  }

  edges <- net %>% distinct()
  edges <- edges %>% inner_join(deseq_results %>%
                                  rownames_to_column("NCBI") %>%
                                  select(from = NCBI,log2FoldChange) %>%
                                  distinct())


  lista_taxones <- data.frame(id = c(edges$from,edges$to)) %>% distinct()

  nodes <- deseq_results %>% rownames_to_column("id") %>% right_join(lista_taxones) %>% left_join(taxtable)

  net2 <- tbl_graph(nodes, edges,node_key = "id")

  begin = min(deseq_results$log2FoldChange)
  end = max(deseq_results$log2FoldChange)

  ggraph(net2, circular =T)  +
    geom_edge_elbow(aes(edge_color=log2FoldChange),edge_width = 10)+
    geom_node_label(aes(label = Name2, fill = log2FoldChange)) + theme_void() +
    scale_fill_viridis(limits=c(begin,end)) +
    scale_edge_color_viridis(limits=c(begin,end))

}


