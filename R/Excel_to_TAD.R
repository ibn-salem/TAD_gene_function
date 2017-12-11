# transforms a excel file with multiple datasheets into a GRangesList object, containing all ranges 
# of every sheet from the file
xltoTAD <- function(exfile, maxsize){
  
  sheet_names <- excel_sheets(exfile)
  gr_list_TAD <- GRangesList()
  
  for(entry in sheet_names){
    
    #reads one sheet of a multisheeted excelfile into a three columned tibble
    tb <- read_xlsx(exfile, sheet = entry, col_names = FALSE)
    
    #Initialisation of counting variables
    i <- 1
    
    #initilaises three vectors to take up chromosome names, starting and endpoints of TADs
    tadStart <- vector("double", nrow(tb))
    tadEnd <- vector("double", nrow(tb))
    tadChrom <- vector("character", nrow(tb))
    
    #iterates over the tibble from first row to second last row
    while(i < nrow(tb)){
      
      #extracts two following chromosome names to compare them
      chromName <- toString(tb[i,1])
      nextChromName <- toString(tb[i+1,1])
      
      # tests wether the starting and end point are on the same chromosome and if so 
      # concatenates each vector
      if (nextChromName != chromName){
        i <- i + 1
      } else {
        #vectors are concatenated with +1 or -1 to exclude the first and last coordinate of
        #the corresponding boundary
        tadStart[i] <- as.numeric(tb[i,3]) + 1
        tadEnd[i] <- as.numeric(tb[i+1,2]) - 1
        tadChrom[i] <- chromName
        i <- i + 1
      }
    }
    #builds a GRanges Object
    gr <- GRanges(tadChrom[tadChrom != ""], IRanges(tadStart[tadStart != "0"], tadEnd[tadEnd != "0"]))
    subgr <- gr[width(gr) <= maxsize]
    grlist <- GRangesList("Celltype" = subgr)
    gr_list_TAD <- c(gr_list_TAD, grlist)
  }
  
  # renames the GRangeslist to its original Celltype names
  names(gr_list_TAD) <- sheet_names
  
  # returns the GRangesList object containing GRanges of all sheets in the entered excel file
  return(gr_list_TAD)
}