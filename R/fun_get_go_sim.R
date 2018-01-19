
geneSimScore <- function(...){
  sim_result <- geneSim(...)
  
  if(is.na(sim_result)){
    return(NA)
  }else{
    return(sim_result$geneSim)
  }
 
}

get_go_sim <- function (g1_id, g2_id, semdata){
  
  map2_dbl(g1_id, g2_id,
    geneSimScore,
    semData = semdata,
    measure = "Resnik",
    combine = "max")

}







