#takes a list of tibbles and a list of their names converts all dataframes in that list 
#into a GRangesList object and returns it
xltoBDY <- function(exfile_sheets, seq_info){
  
  grlist <- GRangesList()
  sheet_names <- c("H1", "IMR90", "GM12878")
  
  for (entry in sheet_names){
    tb <- as_tibble(read_xlsx(exfile_sheets, sheet = entry, col_names = FALSE))
    gr <- GRanges(tb$X__1, 
                IRanges(tb$X__2, tb$X__3),
                seqinfo = seq_info)
    gr <- trim(gr)
    grl <- GRangesList("Celltype" = gr)
    grlist <- c(grlist, grl)
  }
  
  names(grlist) <- sheet_names

  return(grlist)
}
