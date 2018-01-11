jaccard_matrix <- function(grlist){

  #initialise empty vector
  jackoef <- c()

  #get endpoint of for-loop and length of GRangesList
  z <- seq_along(grlist)
  n <- as.numeric(length(grlist))

  #loop over every element in GRangesList
  for (x in z){
    for(y in z){
      #find overlaps between elements of GRangesList
      overlaps <- findOverlaps(grlist[[x]], grlist[[y]])
      
      #calculate the jaccard coefficient to quantify similarity of boundaries

      c <- length(overlaps)
      a_uniq <- sum(!overlapsAny(grlist[[x]], grlist[[y]]))
      b_uniq <- sum(!overlapsAny(grlist[[y]], grlist[[x]]))
      jaccard <- c/ (c + a_uniq + b_uniq )
      jackoef <- c(jackoef, jaccard)
    }
  }

  # construct a matrix with jaccard coeffiecients and heatmap it
  jmat <- matrix(jackoef, n, n)
  celltypes <- names(grlist)
  
  # transfer names of cell types to the according rows and columns in the matrix
  rownames(jmat) <- celltypes
  colnames(jmat) <- celltypes
  

  return(jmat)
} 