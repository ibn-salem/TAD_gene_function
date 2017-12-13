tidy_gosimList <- function(gosimList){
  
  unl_gosim <- unlist(gosimList)
  unl_gosim <- unl_gosim[names(unl_gosim) == "geneSim" | is.na(unl_gosim)]
  unl_gosim <- as.double(unl_gosim)
  
  return(unl_gosim)
}

