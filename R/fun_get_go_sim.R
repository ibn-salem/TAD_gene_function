
get_go_sim <- function (cispair, semdata){
  
   map2_dbl(cispair$g1_id,
            cispair$g2_id,
            geneSim,
            semData = semdata,
            measure = "Resnik",
            combine = "max")
  
}







