#new taxon_polled_I_replicate2 --> taxon_polled_I_replicate2

cI <-  read.table("rawdata/total-spruce-taxa_clade_I.txt",
                  header=TRUE, sep="\t")%>%
  dplyr::select(taxa,clade,value,sample, soil, Cu, waterlevel)

cI <- cI%>%
  mutate(treatment =  paste0(cI$soil,cI$Cu, cI$waterlevel))%>%
  group_by(sample)%>%
  mutate(rel_abund_raw = (value / sum(value))*100)  


cpclade <- Cu_gene%>%
  mutate(sample = as.character(ID),
         cladeI = nosz_cladeI,
         cladeII = nosz_cladeII )%>%
  mutate(sample = str_c("n", sample))%>%
  dplyr::select(sample, cladeI, cladeII)%>%
  mutate(cladeI = as.numeric(cladeI),
         cladeII = as.numeric(cladeII))

str(cI)
cI <- cI%>%
  left_join(cpclade, by = "sample")

cI <- cI%>%
  mutate(cladeIcp_raw = rel_abund_raw*cladeI )







meancI <-  cI %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100,
         N = N)%>%
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund,N, treatment)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
               names_to="level",
               values_to="taxon")

meancI_wide <-  cI %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100)%>%
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund, treatment)


taxon_I_replicate  <- meancI %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 


taxon_poll_I_replicate <- taxon_I_replicate%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")



taxon_polled_I_replicate <- inner_join(taxon_I_replicate, taxon_poll_I_replicate, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(sample, taxon)%>%
  summarise(mean_rel_abund = sum(rel_abund2),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1))





tlist <- meancI_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)


ttaxon_polled_I_replicate <- taxon_polled_I_replicate%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist, by = "sample") 



ttaxon_polled_I_replicate_imputed <- ttaxon_polled_I_replicate %>%
  replace_na(list(" Mesorhizobium" = 0, " Microvirga"  = 0,
                  " Noviherbaspirillum" = 0, " Polaromonas" = 0,
                  " Pseudomonas"   = 0,    " Ramlibacter"  = 0,     " Rhizobium/Agrobacterium_group" = 0, " Skermanella"    = 0,
                  "Other"  = 0,  " unclassified_Acetobacteraceae"  = 0, " Azospirillum"   = 0,  " Neisseria" = 0 ))





microbe <- read_excel("rawdata/microbe_alldata.xlsx")

env <- select(microbe, "ID", "pH",
              "treatment", "soil", "Cu", "waterlevel",
              "replicate", "amoa_Achaea", "amoa_Bacteria", "nosz_cladeI",
              "16s_bacteria", "nosz_cladeII", "Cumass", "NO3", "NH4", "wfps", "flux")%>%
  mutate(ID = as.character(ID))


clade_env_de <- ttaxon_polled_I_replicate2_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "de")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  select(-Cumass, -treatment.y, -replicate, 
         -ID)%>%
  slice(-17:-32)


clade_env1_de  <- clade_env_de %>%
  slice(-17:-32) %>%
  column_to_rownames(var = "sample")%>%
  select(-amoa_Achaea, -amoa_Bacteria,-nosz_cladeI,
         -`16s_bacteria`,
         -nosz_cladeII,-treatment.x, -Cu, -waterlevel,
         -flux)

str(clade_env1_de)


##upland environmental varibels
clade_env_up <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "up")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  # mutate(Cumass = `Cu_mg/kg`)%>%
  select(-Cumass, -treatment.y, -replicate, 
         -ID, -soil,  -amoa_Achaea, -nosz_cladeII)%>%
  slice(-14:-32) 


clade_env1_up <- clade_env_up %>%
  column_to_rownames(var = "sample")%>%
  select( -amoa_Bacteria,-nosz_cladeI,
          -`16s_bacteria`,-treatment.x, -Cu,-waterlevel,
          -flux)




taxon_rel_abundI_genus3 <- meancI %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(rel_abund3 = sum(rel_abund2),
            mean_rel_abund2 = mean(rel_abund3), .groups="drop") %>%
  # mutate(taxon = str_replace(taxon,
  #                            "unclassified_(.*)", "Unclassified *\\1*"),
  #        taxon = str_replace(taxon,
  #                            "^(\\S*)$", "*\\1*"))%>%
  # mutate(taxon = str_replace(taxon,
  #                            "unclassified_(.*)", "Unclassified *\\1*"),
  #        taxon = str_replace(taxon,
  #                            "(.*)_group", "*\\1* group"))%>% #*_group
  # mutate(taxon = str_trim(taxon),
  #        taxon = str_replace_all(taxon, "\\s", ""),  # Remove any whitespace
  #        taxon = if_else(!str_detect(taxon, "^\\*.*\\*$"), paste0("*", taxon, "*"), taxon))%>%
# mutate(taxon = str_replace(taxon, "\\*Unclassified\\*(.*)\\*\\*", "Unclassified *\\1*"))%>%
# mutate(taxon = str_replace(taxon, "\\*\\*(.*)\\*group\\*", "*\\1* group"))%>% 
mutate(replicates = case_when(
  treatment == "dehh" ~ 4,
  treatment == "delh" ~ 4,
  treatment == "dell" ~ 4,
  treatment == "dehl" ~ 4,
  treatment == "uplh" ~ 3,
  treatment == "upll" ~ 4,
  treatment == "uphl" ~ 2,
  treatment == "uphh" ~ 4))%>%
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)





taxon_pollIg3 <- taxon_rel_abundI_genus3%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3,
            mean = mean(mean_rel_abund3),
            .groups = "drop")%>%
  mutate(taxon = str_trim(taxon))

taxon_I_replicate <- taxon_I_replicate%>%
  mutate(taxon = str_trim(taxon))

taxon_polled_I_replicate2 <-
  inner_join(taxon_I_replicate, taxon_pollIg3, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(sample, taxon)%>%
  summarise(mean_rel_abund = sum(rel_abund2),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1))







tlist <- meancI_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)

view(ttaxon_polled_I_replicate2)

ttaxon_polled_I_replicate2 <- taxon_polled_I_replicate2%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist, by = "sample") 
#~asv.tax.t in the example

ttaxon_polled_I_replicate2_imputed <- ttaxon_polled_I_replicate2 %>%
  replace_na(list("Mesorhizobium" = 0, "Microvirga"  = 0,
                  "Noviherbaspirillum" = 0, "Polaromonas" = 0,
                  "Pseudomonas"   = 0,    "Ramlibacter"  = 0,     "Rhizobium/Agrobacterium_group" = 0, "Skermanella" = 0,
                  "Other"  = 0,  "unclassified_Acetobacteraceae"  = 0, "Azospirillum"   = 0,  "Neisseria" = 0 ))









ttaxon_polled_I_replicate2_imputed$sample_for_row_names <- ttaxon_polled_I_replicate2_imputed$sample

# Set the new column as row names
ttaxon_polled_I_replicate2_imputed <- ttaxon_polled_I_replicate2_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate2_imputed)

# Select only numeric columns for metaMDS
numeric_columns <- sapply(ttaxon_polled_I_replicate2_imputed, is.numeric)
ttaxon_numeric2 <- ttaxon_polled_I_replicate2_imputed[, numeric_columns]



##depression soils within both soils
# Set the new column as row names
ttaxon_polled_I_replicate2_imputed_de <- ttaxon_polled_I_replicate2_imputed %>%
  filter(soil == "de")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate2_imputed_de)

# Select only numeric columns for metaMDS
numeric_columns_de <- sapply(ttaxon_polled_I_replicate2_imputed_de, is.numeric)
ttaxon_numeric2_de2 <- ttaxon_polled_I_replicate2_imputed_de[, numeric_columns_de]

##upland soils within both soils
# Set the new column as row names
ttaxon_polled_I_replicate2_imputed_up <- ttaxon_polled_I_replicate2_imputed %>%
  filter(soil == "up")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate2_imputed_up)

# Select only numeric columns for metaMDS
numeric_columns_up <- sapply(ttaxon_polled_I_replicate2_imputed_up, is.numeric)
ttaxon_numeric2_up2 <- ttaxon_polled_I_replicate2_imputed_up[, numeric_columns_up]



# 
# 
# 
# taxon_polled_I_replicate2_up <- inner_join(taxon_I_replicate_up, taxon_poll_I_replicate_up, by = "taxon")%>%
#   mutate(taxon = if_else(pool, "Other", taxon))%>%
#   group_by(sample, taxon)%>%
#   summarise(mean_rel_abund = sum(rel_abund2),
#             mean = min(mean),
#             .groups = "drop")%>%
#   mutate(taxon = factor(taxon),
#          taxon = fct_reorder(taxon, mean, .desc = TRUE),
#          taxon = fct_shift(taxon, n = 1)) %>%
#   #if turn off desc, Firmicutes will be on the bottom
#   mutate(taxon = factor(taxon),
#          taxon = fct_shift(taxon, n = 1))
# 
# 
# 
# ttaxon_polled_I_replicate2_up <- taxon_polled_I_replicate2_up%>%
#   select(-mean)%>%
#   pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
#   right_join( tlist_up, by = "sample")  %>%
#   replace_na(list("Mesorhizobium" = 0, "Microvirga"  = 0,
#                   "Noviherbaspirillum" = 0, "Polaromonas" = 0,
#                   "Pseudomonas"   = 0,    "Ramlibacter"  = 0,     "Rhizobium/Agrobacterium_group" = 0, "Skermanella"    = 0,
#                   "Other"  = 0,  "unclassified_Acetobacteraceae"  = 0, "Azospirillum"   = 0,  "Neisseria" = 0 ))
# 
# 
# 
# 
# # Create a new column for row names
# ttaxon_polled_I_replicate2_up$sample_for_row_names <- ttaxon_polled_I_replicate2_up$sample
# 
# # Set the new column as row names
# ttaxon_polled_I_replicate2_imputed_up <- ttaxon_polled_I_replicate2_up %>%
#   column_to_rownames(var = "sample_for_row_names")
# 
# # Check the structure to identify non-numeric columns
# str(ttaxon_polled_I_replicate2_imputed_up)
# 
# # Select only numeric columns for metaMDS
# numeric_columns_up <- sapply(ttaxon_polled_I_replicate2_imputed_up, is.numeric)
# ttaxon_numeric_up <- ttaxon_polled_I_replicate2_imputed_up[, numeric_columns_up]




###clade II


taxon_rel_abundII_genus2 <- meancII %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(rel_abund3 = sum(rel_abund2),
            mean_rel_abund2 = mean(rel_abund3), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "(.*)_group", "*\\1* group"))%>% #*_group
  mutate(taxon = str_trim(taxon),
         taxon = str_replace_all(taxon, "\\s", ""),  # Remove any whitespace
         taxon = if_else(!str_detect(taxon, "^\\*.*\\*$"), paste0("*", taxon, "*"), taxon))%>%
  mutate(taxon = str_replace(taxon, "\\*Unclassified\\*(.*)\\*\\*", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "\\*\\*(.*)\\*group\\*", "*\\1* group"),
         taxon = str_replace(taxon, "(.*)_(.*)", "\\1 \\2"))%>%
  mutate(replicates = case_when(
    treatment == "dehh" ~ 4,
    treatment == "delh" ~ 4,
    treatment == "dell" ~ 4,
    treatment == "dehl" ~ 4,
    treatment == "uphh" ~ 1))%>%
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)


taxon_rel_abundII_genus2_clean <- meancII %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(rel_abund3 = sum(rel_abund2),
            mean_rel_abund2 = mean(rel_abund3), .groups="drop") %>%
  # mutate(taxon = str_replace(taxon,
  #                            "unclassified_(.*)", "Unclassified *\\1*"),
  #        taxon = str_replace(taxon,
  #                            "(.*)_group", "*\\1* group"))%>% #*_group
  # mutate(taxon = str_trim(taxon),
  #        taxon = str_replace_all(taxon, "\\s", ""),  # Remove any whitespace
  #        taxon = if_else(!str_detect(taxon, "^\\*.*\\*$"), paste0("*", taxon, "*"), taxon))%>%
  # mutate(taxon = str_replace(taxon, "\\*Unclassified\\*(.*)\\*\\*", "Unclassified *\\1*"),
  #        taxon = str_replace(taxon, "\\*\\*(.*)\\*group\\*", "*\\1* group"),
  #        taxon = str_replace(taxon, "(.*)_(.*)", "\\1 \\2"))%>%
  mutate(taxon = str_trim(taxon),
         replicates = case_when(
           treatment == "dehh" ~ 4,
           treatment == "delh" ~ 4,
           treatment == "dell" ~ 4,
           treatment == "dehl" ~ 4,
           treatment == "uphh" ~ 1))%>%
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)

# 
# 
# taxon_pollIIg2 <- taxon_rel_abundII_genus2%>% 
#   group_by(taxon)%>%
#   summarise(pool = max(mean_rel_abund3) <3,
#             mean = mean(mean_rel_abund3),
#             .groups = "drop")
# 
# 
# 

taxon_II_replicate  <- meancII %>%
  filter(soil == "de")%>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100,
         taxon = str_trim(taxon),
         .groups ="drop")



taxon_pollIIg2_2 <- taxon_rel_abundII_genus2_clean%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3,
            mean = mean(mean_rel_abund3),
            .groups = "drop")



taxon_polled_II_replicate2 <- inner_join(taxon_II_replicate, taxon_pollIIg2_2, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(sample, taxon)%>%
  summarise(mean_rel_abund = sum(rel_abund2),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1))



meancII_wide <-  cII %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100)%>%
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund, treatment)




tlist2 <- meancII_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)


ttaxon_polled_II_replicate2 <- taxon_polled_II_replicate2%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist2, by = "sample")


ttaxon_polled_II_replicate_imputed2 <- ttaxon_polled_II_replicate2 %>%
  replace_na(list("Chitinophaga" = 0, "Cytophagaceae_bacterium"  = 0,
                  "Flavisolibacter" = 0, "Opitutaceae_bacterium" = 0,
                  "Pedobacter"   = 0,    "Pontibacter"  = 0,  "Prevotella" = 0, 
                  "Other"  = 0))






ttaxon_polled_II_replicate_imputed2$sample_for_row_names <- ttaxon_polled_II_replicate_imputed2$sample

# Set the new column as row names
ttaxon_polled_II_replicate_imputed2 <- ttaxon_polled_II_replicate_imputed2 %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_II_replicate_imputed2)

# Select only numeric columns for metaMDS
numeric_columnsII2 <- sapply(ttaxon_polled_II_replicate_imputed2, is.numeric)
ttaxon_numeric_II2 <- ttaxon_polled_II_replicate_imputed2[, numeric_columnsII2]


##adonis clade I#####
#with only depressionsoil 
asv.bray1_de2 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_de[,2:(ncol(ttaxon_polled_I_replicate2_imputed_de)-4)], 
                                                 diag = TRUE, upper = TRUE)))

dist_matrix_I_de2 <- as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_de[,2:(ncol(ttaxon_polled_I_replicate2_imputed_de)-4)], 
                                       diag = TRUE, upper = TRUE))


heatmap(dist_matrix_I_de2)

str(dist_matrix_I_de2)

adonis2(asv.bray1_de2 ~ clade_env_de$treatment.x, asv.bray1_de2, permutations = 9000) # p = 0.05688 . treatment
adonis2(asv.bray1_de2 ~ clade_env_de$Cu, asv.bray1_de2, permutations = 9000) # p = 0.01644 * Cu
adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel, asv.bray1_de2, permutations = 9000) # p = 0.3741 waterlevel
adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel*clade_env_de$Cu, asv.bray1_de2, permutations = 9000) # p = 0.01489 *

adonis_I_de2 <- adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel+clade_env_de$Cu, 
                        asv.bray1_de2, permutations = 9000) # p = 0.2266 water level



adonis_I_de_df2 <- data.frame(
  Df = adonis_I_de2$Df,
  SumOfSqs = adonis_I_de2$SumOfSqs,
  R2 = adonis_I_de2$R2,
  F = adonis_I_de2$F,
  `Pr(>F)` = adonis_I_de2$`Pr(>F)`)%>% round(2)

varibles <- c("water", "Cu", "Residual", "total") 

adonis_I_de_df2 <- cbind(varibles, adonis_I_de_df2)

# write_xlsx(adonis_I_de_df2,"~/Downloads/seq/stat/adonis2/adonis_I_de2.xlsx")





#with only upland soils
asv.bray1_up2 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_up[,2:(ncol(ttaxon_polled_I_replicate2_imputed_up)-4)], 
                                                 diag = TRUE, upper = TRUE)))

dist_matrix_I_up2 <- as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_up[,2:(ncol(ttaxon_polled_I_replicate2_imputed_up)-4)], 
                                       diag = TRUE, upper = TRUE))
heatmap(dist_matrix_I_up2)


adonis2(asv.bray1_up2 ~ clade_env_up$Cu, asv.bray1_up2, permutations = 9000) # p = 0.185 Cu
adonis2(asv.bray1_up2 ~ clade_env_up$waterlevel, asv.bray1_up2, permutations = 9000) # p = 0.03266 * Waterlevel
adonis2(asv.bray1_up2 ~ clade_env_up$waterlevel*clade_env_up$Cu, asv.bray1_up2, permutations = 9000) # 0.02511 * Waterlevel
adonis2(asv.bray1_up2 ~ clade_env_up$treatment.x, asv.bray1_up2, permutations = 9000) # p = 0.0711 . treatment

adonis_I_up2 <- adonis2(asv.bray1_up2 ~ clade_env_up$waterlevel+clade_env_up$Cu, asv.bray1_up2, permutations = 9000) # 0.02511 * Waterlevel


adonis_I_up_df2 <- data.frame(
  Df = adonis_I_up2$Df,
  SumOfSqs = adonis_I_up2$SumOfSqs,
  R2 = adonis_I_up2$R2,
  F = adonis_I_up2$F,
  `Pr(>F)` = adonis_I_up2$`Pr(>F)`)%>%round(3)

adonis_I_up_df2 <- cbind(varibles, adonis_I_up_df2)

# write_xlsx(adonis_I_up_df2,"~/Downloads/seq/stat/adonis2/adonis_I_up2.xlsx")



##adonis cladeII####
dist_matrix_II2 <- as.matrix(vegdist(ttaxon_polled_II_replicate_imputed2[,2:(ncol(ttaxon_polled_II_replicate_imputed2)-4)], 
                                     diag = TRUE, upper = TRUE))

heatmap(dist_matrix_II2)

asv.bray1_II2 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_II_replicate_imputed2[,2:(ncol(ttaxon_polled_II_replicate_imputed2)-4)], 
                                                 diag = TRUE, upper = TRUE)))


adonis2(asv.bray1_II2 ~ ttaxon_polled_II_replicate_imputed2$Cu, asv.bray1_II2, permutations = 9000) # p = 0.2663 Cu
adonis2(asv.bray1_II2 ~ ttaxon_polled_II_replicate_imputed2$waterlevel, asv.bray1_II2, permutations = 9000) # p =  0.05155 . waterlevel
adonis2(asv.bray1_II2 ~ ttaxon_polled_II_replicate_imputed2$waterlevel*ttaxon_polled_II_replicate_imputed2$Cu,
        asv.bray1_II2, permutations = 9000) #p = 0.06088  . waterlevel
adonis2(asv.bray1_II2 ~ ttaxon_polled_II_replicate_imputed2$treatment, asv.bray1_II2, permutations = 9000) # P = 0.09566 .

adonis_II2 <- adonis2(asv.bray1_II2 ~ ttaxon_polled_II_replicate_imputed2$waterlevel+ttaxon_polled_II_replicate_imputed2$Cu,
                      asv.bray1_II2, permutations = 9000) #p = 0.06088  . waterlevel


print(adonis_II2)
(adonis_II2$aov.tab) 
as.data.frame(adonis_II2$aov.tab)
adonis_II2_df <- as.data.frame(adonis_II2$aov.tab)
str(adonis_II2)


varibles <- c("water", "Cu", "Residual", "total") 

adonis_II2_df <- data.frame(
  Df = adonis_II2$Df,
  SumOfSqs = adonis_II2$SumOfSqs,
  R2 = adonis_II2$R2,
  F = adonis_II2$F,
  `Pr(>F)` = adonis_II2$`Pr(>F)`)

adonis_II2_df <- cbind(varibles, adonis_II2_df)

# write_xlsx(adonis_II2_df,"~/Downloads/seq/stat/adonis2/adonis_II2.xlsx")


adonis_all2 <- cbind(adonis_I_de_df2, adonis_I_up_df2, adonis_II2_df)

# write_xlsx(adonis_all,"~/Downloads/seq/stat/adonis2/adonis_all2.xlsx")






#envfit flux & pH & NO3 & NH4#####
##de env factors 
clade_env_de <- ttaxon_polled_I_replicate2_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "de")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  select(-Cumass, -treatment.y, -replicate, 
         -ID)%>%slice(-17:-32) 


clade_env1_de_f  <- clade_env_de %>%
  column_to_rownames(var = "sample")%>%
  select(-amoa_Achaea, -amoa_Bacteria,-nosz_cladeI,
         -`16s_bacteria`,
         -nosz_cladeII,-treatment.x,
         -wfps,-soil, -Cu, -waterlevel)



##*de cladeI NMDS####
cladeI_genus_nmds_de2 <- metaMDS(ttaxon_numeric2_de2, k = 2) #distance default is bray

?metaMDS
plot(cladeI_genus_nmds_de2)



envfit_I_de_nmds_f2 <- envfit(cladeI_genus_nmds_de2, clade_env1_de_f, permu = 999, 
                              na.rm = TRUE)



plot(envfit_I_de_nmds_f2, display = "sites")
plot(envfit_I_de_nmds_f2)



envfit_I_de_arrows_nmds_f2 <- as.data.frame(envfit_I_de_nmds_f2$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_de_arrows_nmds_f2 <- envfit_I_de_arrows_nmds_f2%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      TRUE ~ varNames))

# varNames == "Cumass" ~ as.character(c("Cu")),



print(envfit_I_de_arrows_nmds_f2)
# envfit_I_de_r_nmds <- scores(envfit_I_de_nmds$vectors$r)
# envfit_I_de_r_nmds <- sqrt(envfit_I_de_r_nmds)
# envfit_I_de_arrows_nmds$NMDS3  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS4 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1

##significnat test nmds envfit cladeI de####
r_I_de_nmds_env_f2 <- as.matrix(envfit_I_de_nmds_f2$vectors$r) 
p_I_de_nmds_env_f2 <- as.matrix(envfit_I_de_nmds_f2$vectors$pvals)
env_p_I_de_nmds_f2 <- cbind(r_I_de_nmds_env_f2, p_I_de_nmds_env_f2)
colnames(env_p_I_de_nmds_f2) <- c("r2","p_value")
K_I_de_nmds_f2 <- as.data.frame(env_p_I_de_nmds_f2)
K_I_de_nmds_f2 <- K_I_de_nmds_f2%>%
  mutate(p_adjust = p.adjust(K_I_de_nmds_f2$p_value, method = "BH"))
K_I_de_nmds_f2

K_I_de_nmds_f2 <- K_I_de_nmds_f2%>%
  rownames_to_column("de_factors")

# write_xlsx(K_I_de_nmds_f2,"~/Downloads/seq/stat/envfit3/flux/K_I_de_nmds_f2.xlsx")




##envfit cladeI_genus_nmds_de, and plotting######
#sp
cladeI_genus_nmds_points_de_sp2 <- as.data.frame(scores(cladeI_genus_nmds_de2)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_sp12 <- cladeI_genus_nmds_points_de_sp2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_de_si2 <- as.data.frame(scores(cladeI_genus_nmds_de2)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_si12 <- cladeI_genus_nmds_points_de_si2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  





cladeI_genus_nmds_points_de2 <- as.data.frame(cladeI_genus_nmds_de2$points)%>%
  mutate(sample = row.names(.))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))%>%
  slice(-17:-29)

centroid_I_de22 <- cladeI_genus_nmds_points_de2%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )





######de NMDs flux species plotting +envfit #####
# nmds_envfit_plot_I_de 
  cladeI_genus_nmds_points_de2%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  geom_point(data = centroid_I_de22, size = 5, show.legend = TRUE)+
  coord_cartesian(xlim=c(-0.3, 1), ylim=c(-0.6, 1)) +
  geom_segment(data = cladeI_genus_nmds_points_de_sp12,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_de_arrows_nmds_f2,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_de_arrows_nmds_f2,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_de_sp12,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.8,parse = TRUE,
                   size = 3,
                   box.padding = unit(0.5, "lines"),  # Increase padding around label box
                   point.padding = unit(0.5, "lines"), # Increase padding around points
                   nudge_x = 0.15,  # Adjust the nudging for better separation
                   nudge_y = 0.2,
                   segment.size = 0.1) +
  
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  scale_color_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                              depression260Cu60WHC_c1,
                                              depression0Cu90WHC_c1,
                                              depression260Cu90WHC_c1),
                     labels=c("Cu0 WHC60",
                              "Cu260 WHC60",
                              "Cu0 WHC90",
                              "Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  labs(
       x="NMDS1",
       y="NMDS2")+
  # caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
  # envfit(cladeI_genus_nmds_de, clade_env1_de_f, permu = 999, na.rm = TRUE)") +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.75, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")


ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeI_de_nmds_env_species_f2.tiff",  width=9, height=9)
ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeI_de_nmds_env_species_f2.jpeg",  width=9, height=9)






##*upland environmental varibels####
clade_env_up <- ttaxon_polled_I_replicate2_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "up")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  select(-Cumass, -treatment.y, -replicate, 
         -ID,  -amoa_Achaea, -nosz_cladeII)%>%
  slice(-14:-32) 


clade_env1_up_f <- clade_env_up %>%
  column_to_rownames(var = "sample")%>%
  select( -amoa_Bacteria,-nosz_cladeI,
          -`16s_bacteria`,
          -wfps,-soil, -Cu, -waterlevel,-treatment.x)





##upflux cladeI NMDS####
cladeI_genus_nmds_up2 <- metaMDS(ttaxon_numeric2_up2, k = 2) #distance defualt is bray

?metaMDS
plot(cladeI_genus_nmds_up2)



envfit_I_up_nmds_f2 <- envfit(cladeI_genus_nmds_up2, clade_env1_up_f, permu = 999, 
                              na.rm = TRUE)


plot(envfit_I_up_nmds_f2, display = "sites")
plot(envfit_I_up_nmds_f2)





envfit_I_up_arrows_nmds_f2 <- as.data.frame(envfit_I_up_nmds_f2$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_up_arrows_nmds_f2 <- envfit_I_up_arrows_nmds_f2%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_I_up_arrows_nmds2)
envfit_I_up_r_nmds <- scores(envfit_I_up_nmds_f2$vectors$r)
envfit_I_up_r_nmds <- sqrt(envfit_I_up_r_nmds)
# envfit_I_up_arrows_nmds$NMDS3  <- envfit_I_up_arrows_nmds$NMDS1  / as.numeric(envfit_I_up_r_nmds) * 1
# envfit_I_up_arrows_nmds$NMDS4 <- envfit_I_up_arrows_nmds$NMDS2 / as.numeric(envfit_I_up_r_nmds) * 1

##significnat test nmds envfit cladeI up####
r_I_up_nmds_env_f2 <- as.matrix(envfit_I_up_nmds_f2$vectors$r) 
p_I_up_nmds_env_f2 <- as.matrix(envfit_I_up_nmds_f2$vectors$pvals)
env_p_I_up_nmds_f2 <- cbind(r_I_up_nmds_env_f2, p_I_up_nmds_env_f2)
colnames(env_p_I_up_nmds_f2) <- c("r2","p_value")
K_I_up_nmds_f2 <- as.data.frame(env_p_I_up_nmds_f2)%>%
  mutate(p_adjust = p.adjust(K_I_up_nmds_f2$p_value, method = "BH"))
K_I_up_nmds_f2

K_I_up_nmds_f2 <- K_I_up_nmds_f2%>% 
  rownames_to_column(var = "up_facotrs")

# write_xlsx(K_I_up_nmds_f2,"~/Downloads/seq/stat/envfit3/flux/K_I_up_nmds_f2.xlsx")



##envfit cladeI_genus_nmds_up, and plotting######
#sp
cladeI_genus_nmds_points_up_sp2 <- as.data.frame(scores(cladeI_genus_nmds_up2)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_sp12 <- cladeI_genus_nmds_points_up_sp2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_up_si2 <- as.data.frame(scores(cladeI_genus_nmds_up2)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_si12 <- cladeI_genus_nmds_points_up_si2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  



cladeI_genus_nmds_points_up2 <- as.data.frame(cladeI_genus_nmds_up2$points)%>%
  mutate(sample = row.names(.))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "upll",
                                "uphl",
                                "uplh",
                                "uphh")))%>%
  slice(-14:-29)


centroid_I_up22 <- cladeI_genus_nmds_points_up2%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )


######up NMDs species plotting +envfit #####
# cladeI_genus_nmds_points_up_sp%>%
# nmds_envfit_plot_I_up <- 
  cladeI_genus_nmds_points_up2%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  coord_cartesian(xlim=c(-1, 2.3), ylim=c(-1.5, 1)) +
  geom_point(data = centroid_I_up22, size = 5, show.legend = TRUE)+
  geom_segment(data = cladeI_genus_nmds_points_up_sp12,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_up_arrows_nmds_f2,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_up_arrows_nmds_f2,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_up_sp12,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE,
                   size = 3.5,
                   box.padding = unit(0.5, "lines"),  # Increase padding around label box
                   point.padding = unit(0.5, "lines"), # Increase padding around points
                   nudge_x = 0,  # Adjust the nudging for better separation
                   nudge_y = 0.1,
                   segment.size = 0.1) +
  scale_color_manual(name = NULL, values  = c(upland0Cu60WHC_c1,
                                              upland260Cu60WHC_c1,
                                              upland0Cu90WHC_c1,
                                              upland260Cu90WHC_c1),
                     labels=c("Cu0 WHC60",
                              "Cu260 WHC60",
                              "Cu0 WHC90",
                              "Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(upland0Cu60WHC_c1,
                                             upland260Cu60WHC_c1,
                                             upland0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  labs(
       x="NMDS1",
       y="NMDS2")+
  # caption = "Upland soils, using upland soil data relative abundance for metaMDS(k = 2),
  # envfit(cladeI_genus_nmds_up, clade_env1_up_f, permu = 999, na.rm = TRUE)") +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")


ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeI_up_nmds_env_species_f2.tiff",  width=9, height=9)
ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeI_up_nmds_env_species_f2.jpeg",  width=9, height=9)




##*with flux clade II#####


cladeII_genus_nmds2 <- metaMDS(ttaxon_numeric_II2, k = 2)

envfit_II_de_nmds_f2 <- envfit(cladeII_genus_nmds2, clade_env1_de_f, permu = 999, 
                               na.rm = TRUE)

r_II_de_nmds_env_f2 <- as.matrix(envfit_II_de_nmds_f2$vectors$r) 
p_II_de_nmds_env_f2 <- as.matrix(envfit_II_de_nmds_f2$vectors$pvals)
env_p_II_de_nmds_f2 <- cbind(r_II_de_nmds_env_f2,p_II_de_nmds_env_f2)
colnames(env_p_II_de_nmds_f2) <- c("r2","p_value")
K_II_de_nmds_f2 <- as.data.frame(env_p_II_de_nmds_f2)%>%
  mutate(p_adjust = p.adjust(K_II_de_nmds_f2$p_value, method = "BH"))
K_II_de_nmds_f2


K_II_de_nmds_f2 <- K_II_de_nmds_f2%>% 
  rownames_to_column(var = "II_facotrs")

write_xlsx(K_II_de_nmds_f2,"~/Downloads/seq/stat/envfit3/flux/K_II_de_nmds_f2.xlsx")



K_I_II_all <- cbind(K_I_de_nmds_f2, K_I_up_nmds_f2,K_II_de_nmds_f2)

write_xlsx(K_I_II_all,"~/Downloads/seq/stat/envfit3/flux/K_I_II_all.xlsx")


envfit_II_de_arrows_nmds_f2 <- as.data.frame(envfit_II_de_nmds_f2$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_II_de_arrows_nmds_f2 <- envfit_II_de_arrows_nmds_f2%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_II_de_arrows_nmds_f2)
# envfit_II_de_r_nmds_f2 <- scores(envfit_II_de_nmds_f2$vectors$r)
# envfit_II_de_r_nmds_f2 <- sqrt(envfit_II_de_r_nmds_f2)
# envfit_I_de_arrows_nmds$NMDS1  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS2 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1





##envfit cladeII_genus_nmds_de, and plotting######
#sp
cladeII_genus_nmds_points_de_sp2 <- as.data.frame(scores(cladeII_genus_nmds2)$species)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_sp12 <- cladeII_genus_nmds_points_de_sp2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_(.*)", "\\1 \\2"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeII_genus_nmds_points_de_si2 <- as.data.frame(scores(cladeII_genus_nmds2)$site)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_si12 <- cladeII_genus_nmds_points_de_si2 %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  





cladeII_genus_nmds_points2 <- as.data.frame(cladeII_genus_nmds2$points)%>%
  mutate(sample = row.names(cladeII_genus_nmds2$points))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist2, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))%>%
  slice(-14:-29)


centroid_II_22 <- cladeII_genus_nmds_points2%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )


######de NMDs species plotting +envfit #####
# cladeII_genus_nmds_points_de_sp%>%
cladeII_genus_nmds_points2%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  geom_point(data = centroid_II_22, size = 5, show.legend = TRUE)+
  coord_cartesian(xlim=c(-1, 1.5), ylim=c(-1, 1)) +
  geom_segment(data = cladeII_genus_nmds_points_de_sp12,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_II_de_arrows_nmds_f2,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_II_de_arrows_nmds_f2,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeII_genus_nmds_points_de_sp12,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE,
                   size = 3.5,
                   box.padding = unit(0.5, "lines"),  # Increase padding around label box
                   point.padding = unit(0.5, "lines"), # Increase padding around points
                   nudge_x = 0.1,  # Adjust the nudging for better separation
                   nudge_y = 0.1,
                   segment.size = 0.1) +
  
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  scale_color_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                              depression260Cu60WHC_c1,
                                              depression0Cu90WHC_c1,
                                              depression260Cu90WHC_c1),
                     labels=c("Cu0 WHC60",
                              "Cu260 WHC60",
                              "Cu0 WHC90",
                              "Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(
       x="NMDS1",
       y="NMDS2"
  ) + 
  # caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
  # envfit(cladeII_genus_nmds_de, clade_env1_de_f, permu = 999, na.rm = TRUE)"
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")


ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeII_de_nmds_env_species_flux2.tiff",  width=9, height=9)
ggsave("~/Downloads/seq/envfitplot_dot/final_only_flux/cladeII_de_nmds_env_species_flux2.jpeg",  width=9, height=9)








taxon_rel_abundI_genus2 <- meancI %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(rel_abund3 = sum(rel_abund2),
            mean_rel_abund2 = mean(rel_abund3), .groups="drop") %>%
  # mutate(taxon = str_replace(taxon,
  #                            "unclassified_(.*)", "Unclassified *\\1*"),
  #        taxon = str_replace(taxon,
  #                            "^(\\S*)$", "*\\1*"))%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "(.*)_group", "*\\1* group"))%>% #*_group
  mutate(taxon = str_trim(taxon),
         taxon = str_replace_all(taxon, "\\s", ""),  # Remove any whitespace
         taxon = if_else(!str_detect(taxon, "^\\*.*\\*$"), paste0("*", taxon, "*"), taxon))%>%
  mutate(taxon = str_replace(taxon, "\\*Unclassified\\*(.*)\\*\\*", "Unclassified *\\1*"))%>%
  mutate(taxon = str_replace(taxon, "\\*\\*(.*)\\*group\\*", "*\\1* group"))%>% 
  mutate(replicates = case_when(
    treatment == "dehh" ~ 4,
    treatment == "delh" ~ 4,
    treatment == "dell" ~ 4,
    treatment == "dehl" ~ 4,
    treatment == "uplh" ~ 3,
    treatment == "upll" ~ 4,
    treatment == "uphl" ~ 2,
    treatment == "uphh" ~ 4))%>%
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))



taxon_pollIg2 <- taxon_rel_abundI_genus2%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3,
            mean = mean(mean_rel_abund3),
            .groups = "drop")
print(taxon_pollIg2, n = 64)
library(RColorBrewer)

colors_set3 <- brewer.pal(12, "Set3")
colors_paired <- brewer.pal(1, "Paired")  # Adding one color from another palette

# Combine the colors to have a total of 13
all_colors <- c(colors_set3, colors_paired)
all_colors <- c("#8DD3C7", "#FFFFB3",  "#1F78B4",  "#A6CEE3", "#FB8072", "#7570B3",
                "#FDB462" ,"#B3DE69" ,"#FCCDE5", "#BC80BD" , "#FFED6F", "#E5C494","gray","#FED9A6")

# "#CCEBC5", "#BEBADA",


inner_join(taxon_rel_abundI_genus2, taxon_pollIg2, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c("upll",
                                "uphl",
                                "uplh",
                                "uphh",
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))%>%
  group_by(treatment, taxon)%>%
  summarise(mean_rel_abund = sum(mean_rel_abund3),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon,mean, .desc = TRUE), 
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1),
         soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"),
         treatment1 = str_replace(treatment,"^\\w*(\\w)(\\w)", "\\1\\2")) %>%
  ggplot(aes(x=treatment, y=mean_rel_abund, fill=taxon)) +
  # facet_wrap(.~soil)+
  geom_col() +
  scale_fill_manual(name=NULL,
                    breaks=c("*Afipia*", "*Azospirillum*","*Bradyrhizobium*","*Hyphomicrobium*",
                             "*Mesorhizobium*", "*Microvirga*", "*Noviherbaspirillum*",
                             "*Pseudomonas*","*Ramlibacter*","*Rhizobium/Agrobacterium* group",
                             "*Skermanella*", "Unclassified *Acetobacteraceae*",
                             "Other"),
                    values = all_colors) +
  #the last color "gray" --> "Other"
  # scale_fill_manual(name = NULL, 
  #                   
  #                   values = all_colors) +
  scale_x_discrete(breaks=c("dehh",
                            "dehl",
                            "delh",
                            "dell",
                            "uphh",
                            "uphl",
                            "uplh",
                            "upll"),
                   labels = c("depression<br>Cu260<br>water90",
                              "depression<br>Cu260<br>water60",
                              "depression<br>Cu0<br>water90",
                              "depression<br>Cu0<br>water60",
                              "upland<br>Cu260<br>water90",
                              "upland<br>Cu260<br>water60",
                              "upland<br>Cu0<br>water90",
                              "upland<br>Cu0<br>water60"))+
  # labels=c("Healthy",
  #          "Diarrhea,<br>*C. difficile*<br>negative",
  #          "Diarrhea,<br>*C. difficile*<br>positive")) +
  # scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        axis.text.y = element_markdown(size = 12),
        legend.text = element_markdown(size = 12),
        legend.key.size = unit(28, "pt"),
        axis.title.y  = element_markdown(size = 12),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
# annotate("text", x = 2.5, y = -3, label = "Upland soil", size = 5, angle = 0, hjust = 0.5) +
# annotate("text", x = 6.5, y = -3, label = "Depression soil", size = 5, angle = 0, hjust = 0.5) 

