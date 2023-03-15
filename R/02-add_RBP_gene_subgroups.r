# Aim: to include gene category information in the characterisation of gene expression variation across tissues

# function to wrangle gene subgroup xlsx file
wrangle.gene.excel = function(file){
  
  library(magrittr)
  library(janitor)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(tidyverse)
  
  # file="./RBPs_subgroups.xlsx"
  
  # get names of subgroup cols
  cols = readxl::read_excel(file) %>% 
    janitor::clean_names() %>% 
    dplyr::select(-all_of(c("name", "id"))) %>% 
    colnames() 
  
  # wrangle .xlsx file to usable format
  rbp.subgroups = readxl::read_excel(file) %>% 
    janitor::clean_names() %>% 
    # replace 0|1 notation with string denoting gene subgroup
    dplyr::mutate(
      splicing_regulation = case_when(
        splicing_regulation==1 ~ "splicing_regulation"),
      spliceosome = case_when(
        spliceosome==1 ~ "spliceosome"),
      exon_junction_complex = case_when(
        exon_junction_complex==1 ~ "exon_junction_complex")
    ) %>%
    # use data.table to force the 3 subcategory cols into a list
    data.table::as.data.table() %>% 
    .[ , subgroup.listcol := lapply(transpose(.SD), c), .SDcols = cols] %>% 
    # then force back to dataframe
    as.data.frame() %>% 
    tibble::tibble() %>% 
    # remove sep. subgroup cols
    dplyr::select(-all_of(cols)) %>% 
    # unnest to get a sep row for each entry in list in subgroup.listcol
    unnest(cols=c(subgroup.listcol)) %>% 
    # unnest to transform col back to a flat char col 
    unnest(cols=c(subgroup.listcol)) %>% 
    tidyr::drop_na() %>% 
    dplyr::rename(gene.subgroup=subgroup.listcol)
  
  return(rbp.subgroups)
}

# function to plot per-subgroup heatmaps displaying average TPM
wrangle.data.gen.heatmap.avgTPM = function(data.path, rbp.subgroup.xlsx="./RBPs_subgroups.xlsx", txt.size=10){
  
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
        tidyr::pivot_longer(cols=2:ncol(.), names_to="sample.id", values_to="tpm") %>% 
        # generate participant.id 
        dplyr::mutate(participant.id = stringr::str_match_all(string=sample.id, pattern="(GTEX-.*?-.*?-SM)-")[[1]][,2]) %>% 
        # assign file names (tissue names) to column
        dplyr::mutate(tissue.name = fs::path_file(x)) %>% 
        return(.)
    }) %>% 
    dplyr::bind_rows() %>% 
    relocate(participant.id, sample.id, tissue.name, gene.id, tpm)
  
  # get RBP subgroup information
  rbp.subgroups = wrangle.gene.excel(rbp.subgroup.xlsx)
  
  # map RBP subgroup information to expression df
  tissue.dat = tissue.dat %>% 
    dplyr::filter(gene.id %in% (rbp.subgroups$id %>% unique())) %>% 
    dplyr::left_join(x=., y=rbp.subgroups, by=c("gene.id"="id"))
  
  # gene X tissue heatmap:
  p=tissue.dat %>% 
    dplyr::group_by(tissue.name, gene.id) %>% 
    dplyr::mutate(median.tpm = median(tpm, na.rm=T)) %>% 
    dplyr::mutate(median.tpm = as.numeric(median.tpm)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(gene.subgroup) %>% 
    dplyr::group_split() %>% 
    lapply(X=., FUN=function(x){
      return(ggplot(data=x, aes(x=tissue.name, y=reorder(name, median.tpm), fill=median.tpm)) + 
               theme_bw(base_size=txt.size) + 
               geom_tile() + 
               scale_fill_gradient2(low="white", high="#E64B35B2", limits=c(0,650)) + 
               theme(axis.text.x = element_text(angle=45, hjust=1), 
                     legend.position = "top") +
               labs(fill="Median TPM") + 
               facet_grid(~gene.subgroup) + 
               xlab("") + 
               ylab(""))
    })
  
  fig = ( p[[1]] | p[[2]] | p[[3]] ) + plot_layout(guides = 'collect')
  
  return(fig)
  
}

# function to plot per-subgroup heatmaps displaying log-fold change TPM
wrangle.data.gen.heatmap.lfc = function(data.path, rbp.subgroup.xlsx="./RBPs_subgroups.xlsx", txt.size=10){
  
  # data.path="./data/RBP/"
  # rbp.subgroup.xlsx="./RBPs_subgroups.xlsx"
  # txt.size=6

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
        tidyr::pivot_longer(cols=2:ncol(.), names_to="sample.id", values_to="tpm") %>%
        # generate participant.id
        dplyr::mutate(participant.id = stringr::str_match_all(string=sample.id, pattern="(GTEX-.*?-.*?-SM)-")[[1]][,2]) %>%
        # assign file names (tissue names) to column
        dplyr::mutate(tissue.name = fs::path_file(x)) %>%
        return(.)
    }) %>%
    dplyr::bind_rows() %>%
    relocate(participant.id, sample.id, tissue.name, gene.id, tpm)
  
  # get RBP subgroup information
  rbp.subgroups = wrangle.gene.excel(rbp.subgroup.xlsx)
  
  # map RBP subgroup information to expression df
  tissue.dat = tissue.dat %>%
    dplyr::filter(gene.id %in% (rbp.subgroups$id %>% unique())) %>%
    dplyr::left_join(x=., y=rbp.subgroups, by=c("gene.id"="id"))
  
  # gene X tissue heatmap:
  p=tissue.dat %>%
    dplyr::group_by(gene.subgroup) %>%
    dplyr::group_split() %>%
    lapply(X=., FUN=function(x){
      
      # calculate fold change for each subgroup
      x %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(tpm>0) %>%
        dplyr::group_by(tissue.name, gene.id) %>% 
        dplyr::mutate(med.tissue.gene = median(tpm, na.rm=T)) %>% # median across samples, but per tissue-gene pair
        dplyr::ungroup() %>% 
        dplyr::select(-participant.id, -sample.id, -tpm) %>% 
        dplyr::distinct() %>% 
        dplyr::group_by(gene.id) %>% 
        dplyr::mutate(med.med.tissue = median(med.tissue.gene, na.rm=T)) %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(med.med.tissue) %>% 
        dplyr::mutate(label = case_when(
          med.tissue.gene==med.med.tissue ~ "*",
          TRUE ~ ""
        )) %>% 
        # add fold change - using middle tissue as denom. 
        dplyr::mutate(tpm.logfc = log2((med.tissue.gene)/(med.med.tissue))) %>%
        
        ggplot(data=., aes(x=reorder(tissue.name, med.tissue.gene), y=reorder(name, med.med.tissue), fill=tpm.logfc)) +
        theme_bw(base_size=txt.size) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid='white', high="#E64B35B2", limits=c(-10,10)) +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              legend.position = "top") +
        labs(fill="Log-fold change\nwrt med tissue") +
        facet_grid(~gene.subgroup) +
        xlab("") +
        ylab("") %>%
        return(.)
      
    })
  
  fig = ( p[[1]] | p[[2]] | p[[3]] ) + plot_layout(guides = 'collect')
  
  return(fig)
  
}
