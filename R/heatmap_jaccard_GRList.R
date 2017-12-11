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
      a <- length(grlist[[x]])
      b <- length(grlist[[y]])
      c <- length(overlaps)
      jaccard <- c/ sum(a,b)
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