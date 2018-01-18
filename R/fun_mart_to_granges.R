mart_to_granges <- function(mar, seqinfo){
  
  #filters the df to exclude all chromosomes without proper names
  mar <- mar %>% 
    filter(chromosome_name %in% paste0(c(as.character(1:22), "X"))) %>% 
    mutate("transcript_size" = transcript_end - transcript_start) %>% 
    arrange(ensembl_gene_id, desc(transcript_size)) %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
    arrange(chromosome_name, transcription_start_site)
  
  #gets all gene pairs from genes with GO Term with a distance <= 1000000 bp
  granges <- GRanges(seqnames = paste0("chr" ,mar$chromosome_name),
                     IRanges(start = mar$transcription_start_site,
                             end = mar$transcription_start_site),
                     seqinfo = seqinfo,
                     "gene_id" = mar$ensembl_gene_id,
                     "go_term" = mar$go_id,
                     "go_linkage_type" = mar$go_linkage_type,
                     "NCBI_ID" = paste0("hsa:", mar$entrezgene)
                      )
    
  return(granges)
}