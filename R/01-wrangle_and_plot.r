wrangle.data.gen.heatmap = function(data.path, gene.list, fig.file.name, txt.size=10){
  
  # test
  # gene.list = read_delim("./NMD.txt", delim="\t", col_names=T)$Name %>%
  #   ensembldb::select(EnsDb.Hsapiens.v79, keys=., keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
  # data.path="./data/NMD"
  # fig.file="NMD_plot"

  library(tidyverse)
  library(parallel)
  library(stringr)
  library(ensembldb)
  library(EnsDb.Hsapiens.v79)
  library(ggsci)
  library(RColorBrewer)
  library(patchwork)
  library(magrittr)
  
  # get all tissues, and bind into one df
  tissue.dat = list.dirs(path=data.path, full.names=T, recursive=F) %>% 
    parallel::mclapply(X=., mc.cores=10, FUN=function(x){
      vroom::vroom(file=paste0(x,"/tpm_ensembl105.csv")) %>%
        dplyr::rename(gene.id=`...1`) %>% 
        # dplyr::filter(rowSums(across(where(is.numeric))) > 0) %>% # only keep genes with sum(tpm)>0
        # dplyr::select_if(~ !is.numeric(.) || sum(.) > 0) %>% # only keep samples with sum(tpm)>0
        tidyr::pivot_longer(cols=2:ncol(.), names_to="sample.id", values_to="tpm") %>% 
        # generate participant.id 
        dplyr::mutate(participant.id = stringr::str_match_all(string=sample.id, pattern="(GTEX-.*?-.*?-SM)-")[[1]][,2]) %>% 
        # assign file names (tissue names) to column
        dplyr::mutate(tissue.name = fs::path_file(x)) %>% 
        return(.)
    }) %>% 
    dplyr::bind_rows() %>% 
    relocate(participant.id, sample.id, tissue.name, gene.id, tpm) %>% 
    # filter for gene list
    dplyr::filter(gene.id %in% gene.list$GENEID) %>% 
    dplyr::filter()
  
  # gene x tissue heatmap - avg tpm 
  p=tissue.dat %>% 
    dplyr::group_by(tissue.name, gene.id) %>% 
    dplyr::summarise(median.tpm = median(tpm, na.rm=T)) %>% 
    dplyr::mutate(median.tpm=as.numeric(median.tpm)) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(x=., y=gene.list, by=c("gene.id"="GENEID")) %>% 
    
    ggplot(., aes(x=reorder(tissue.name, median.tpm), y=reorder(SYMBOL, median.tpm), fill=median.tpm)) + 
    theme_bw(base_size=txt.size) + 
    geom_tile() + 
    scale_fill_gradient2(low="white", high="#E64B35B2") + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(fill="Median TPM") + 
    xlab("") + 
    ylab("")
  
  # stand dev boxplot
  q=tissue.dat %>% 
    dplyr::group_by(gene.id, participant.id) %>% 
    dplyr::summarise(sd = sd(tpm, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(x=., y=gene.list, by=c("gene.id"="GENEID")) %>% 
    
    ggplot(data=., aes(x=reorder(SYMBOL, sd), y=sd)) + 
    theme_bw(base_size=txt.size) + 
    geom_boxplot(outlier.shape = NA, alpha=0.8) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    xlab("") + 
    ylab("cross-tissue sd") + 
    coord_flip()
  
  # tissue variance boxplot
  r=tissue.dat %>% 
    dplyr::group_by(gene.id, participant.id) %>% 
    dplyr::summarise(variance = var(tpm, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(x=., y=gene.list, by=c("gene.id"="GENEID")) %>% 
    
    ggplot(data=., aes(x=reorder(SYMBOL, variance), y=variance)) + 
    theme_bw(base_size=txt.size) + 
    geom_boxplot(outlier.shape = NA, alpha=0.6) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    xlab("") + 
    ylab("cross-tissue variance") + 
    coord_flip()
  
  # lfc wrt the median tissue
  s=tissue.dat %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(tpm>0) %>%
    dplyr::group_by(tissue.name, gene.id) %>% 
    dplyr::mutate(med.tissue.gene = median(tpm, na.rm=T)) %>% # median across samples, but per tissue-genepair
    dplyr::ungroup() %>% 
    dplyr::select(-participant.id, -sample.id, -tpm) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(gene.id) %>% 
    dplyr::mutate(med.med.tissue = median(med.tissue.gene, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(med.tissue.gene, med.med.tissue) %>% 
    dplyr::mutate(label = case_when(
      med.tissue.gene==med.med.tissue ~ "*",
      TRUE ~ ""
    )) %>% 
    # add fold change - using middle tissue as denom. 
    dplyr::mutate(tpm.logfc = log2((med.tissue.gene)/(med.med.tissue))) %>% 
    dplyr::left_join(x=., y=gene.list, by=c("gene.id"="GENEID")) %>% 
    
    ggplot(., aes(x=tissue.name, y=reorder(SYMBOL, med.med.tissue), fill=tpm.logfc)) + 
    theme_bw(base_size=txt.size) + 
    geom_tile() + 
    scale_fill_gradient2(low="blue", mid='white', high="#E64B35B2") + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(fill="Log-fold change\nwrt med tissue") + 
    # geom_text(mapping=aes(label=label), size=2) +
    xlab("") + 
    ylab("")
  
  plot = ( (p/q) | (r/s) ) + plot_annotation(tag_levels = 'A')
  
  # ggsave(
  #   filename=paste0("./img/", fig.file.name, ".png"),
  #   plot = last_plot(),
  #   device ="png",
  #   scale = 1,
  #   width = 20,
  #   height = 15,
  #   units = c("in"),
  #   dpi = 300
  # )
  
  return(plot)
  
}