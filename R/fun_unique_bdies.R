#creates tibble out of granges list with all boundaries and counts how often every
#start coordinate occurs.

#TODO not sure whether this is sufficient to only use the starting point of a boundary,
#since it is 40001 bp wide
bdies_tibble <- as_tibble(schmitt_BDY)
all_bdies_gr <- unlist(schmitt_BDY)

hit_tibble <- plyr::count(bdies_tibble, "start", wt_var = NULL)
grs <- left_join(bdies_tibble, hit_tibble, by = "start")

#counts the number of boundaries that are uniqe for a celltype and plots it in a barplot
unique_bdies <- grs %>% 
  filter(freq == 1) %>% 
  group_by(group_name) %>% 
  summarise(n = n())

ggplot(unique_bdies, aes( x = group_name, y = n)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Celltype", y = "Nr of uniqe Boundaries")




olList <- map(as.list(schmitt_BDY), ~overlapsAny(all_bdies_gr, .x))
olDF <- as.tibble(olList)

bdies_tibble <- bind_cols(bdies_tibble, olDF)

bdies_tibble <- bdies_tibble %>% 
  select(-c(width, strand, group, group_name))

bdies_tibble <- bdies_tibble %>% 
  mutate(bdy_nr = seq_along(1:nrow(bdies_tibble)))%>% 
  gather(key = "celltype", value = "in_cell", GM12878:SX)

bdies_tibble %>% 
  group_by(bdy_nr) %>% 
  count("in_cell" == TRUE)


  
#TODO create a Grangesobject with resized boundaries (start-500kb, end+500kb) to find
#Genes that overlap with these new ranges


bdies_tibble %>% 
  group_by(seqnames) %>% 
  summarise(n = n()) %>% 
  ggplot(aes( x = seqnames, y = n)) +
  geom_bar( stat = "identity") +
  theme_minimal()+
  labs(x = " ", y = "Nr of boundaries")







