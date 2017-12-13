get_go_sim <- function (cispair, semdata){
  
  map2(cispair$g1_id,
       cispair$g2_id,
       geneSim,
       semData = semdata,
       measure = "Resnik",
       combine = "max")
  
}