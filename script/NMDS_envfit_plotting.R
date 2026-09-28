##gene composition####
taxon_rel_abundI_genus2%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))%>%
  group_by( soil, taxon) %>%
  summarise(N = n(), 
            sum_abund =  sum(mean_rel_abund3),
            average_abund = sum_abund/N,
            .groups = 'drop')%>%
  arrange(desc(sum_abund))



taxon_rel_abundII_genus2%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))%>%
  group_by( soil, taxon) %>%
  summarise(N = n(), 
            sum_abund =  sum(mean_rel_abund3),
            average_abund = sum_abund/N,
            .groups = 'drop')%>%
  arrange(desc(sum_abund))





###pool rel_abund <3 for de and up sepreately#####
tlist <- meancI_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)

tlist_de <- tlist%>%
  filter(soil == "de")

tlist_up <- tlist%>%
  filter(soil == "up")


##data only for depression soils
# taxon_I_replicate <- taxon_I_replicate%>%
#   mutate(taxon = str_trim(taxon))
# 
# taxon_I_replicate_de <- taxon_I_replicate%>%
#   filter(soil == "de")
# 
# taxon_poll_I_replicate_de <- taxon_I_replicate_de%>%
#   group_by(taxon)%>%
#   summarise(pool = max(rel_abund2) <3,
#             mean = mean(rel_abund2),
#             .groups = "drop")
# 
# 
# 
# taxon_polled_I_replicate_de <- inner_join(taxon_I_replicate_de, taxon_poll_I_replicate_de, by = "taxon")%>%
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
# ttaxon_polled_I_replicate_de <- taxon_polled_I_replicate_de%>%
#   select(-mean)%>%
#   pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
#   right_join( tlist_de, by = "sample") 
# 
# 
# 
# ######depression soil by treatment cladeI_genus_nmds_points_de with rel_abundance
# # Create a new column for row names
# ttaxon_polled_I_replicate_de$sample_for_row_names <- ttaxon_polled_I_replicate_de$sample
# 
# # Set the new column as row names
# ttaxon_polled_I_replicate_imputed_de <- ttaxon_polled_I_replicate_de %>%
#   column_to_rownames(var = "sample_for_row_names")
# 
# # Check the structure to identify non-numeric columns
# str(ttaxon_polled_I_replicate_imputed_de)
# 
# # Select only numeric columns for metaMDS
# numeric_columns_de <- sapply(ttaxon_polled_I_replicate_imputed_de, is.numeric)
# ttaxon_numeric_de <- ttaxon_polled_I_replicate_imputed_de[, numeric_columns_de]

##data used before
ttaxon_numeric_de2
# 
# ##data only for up soils
# taxon_I_replicate_up <- taxon_I_replicate%>%
#   filter(soil == "up")
# 
# taxon_poll_I_replicate_up <- taxon_I_replicate_up%>%
#   group_by(taxon)%>%
#   summarise(pool = max(rel_abund2) <3,
#             mean = mean(rel_abund2),
#             .groups = "drop")



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



 meancI %>%
  filter(level=="genus") %>%
  group_by(treatment, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(treatment, taxon) %>%
  summarize(mean_rel_abund = mean(rel_abund), .groups="drop") %>%
  filter(!is.na(taxon)) %>%
  mutate(taxon = str_trim(taxon),
         taxon = str_replace_all(taxon, "\\s", ""),  # Remove any whitespace
         taxon = if_else(!str_detect(taxon, "^\\*.*\\*$"), paste0("*", taxon, "*"), taxon))



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
   mutate(mean_rel_abund3 = mean_rel_abund2/replicates)
 
 
 
 taxon_pollIg2 <- taxon_rel_abundI_genus2%>% 
   group_by(taxon)%>%
   summarise(pool = max(mean_rel_abund3) <3,
             mean = mean(mean_rel_abund3),
             .groups = "drop")
 
 
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
         taxon = fct_shift(taxon, n = 1)) 


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

view(ttaxon_polled_I_replicate)
ttaxon_polled_I_replicate <- taxon_polled_I_replicate%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist, by = "sample") 
#~asv.tax.t in the example

ttaxon_polled_I_replicate_imputed <- ttaxon_polled_I_replicate %>%
  replace_na(list(" Mesorhizobium" = 0, " Microvirga"  = 0,
                  " Noviherbaspirillum" = 0, " Polaromonas" = 0,
                  " Pseudomonas"   = 0,    " Ramlibacter"  = 0,     " Rhizobium/Agrobacterium_group" = 0, " Skermanella"    = 0,
                  "Other"  = 0,  " unclassified_Acetobacteraceae"  = 0, " Azospirillum"   = 0,  " Neisseria" = 0 ))






ttaxon_polled_I_replicate_imputed$sample_for_row_names <- ttaxon_polled_I_replicate_imputed$sample

# Set the new column as row names
ttaxon_polled_I_replicate_imputed <- ttaxon_polled_I_replicate_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed)

# Select only numeric columns for metaMDS
numeric_columns <- sapply(ttaxon_polled_I_replicate_imputed, is.numeric)
ttaxon_numeric <- ttaxon_polled_I_replicate_imputed[, numeric_columns]

ttaxon_numeric

##depression soils within both soils
# Set the new column as row names
ttaxon_polled_I_replicate_imputed_de <- ttaxon_polled_I_replicate_imputed %>%
  filter(soil == "de")%>%
  select(-` Neisseria`)%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed_de)

# Select only numeric columns for metaMDS
numeric_columns_de <- sapply(ttaxon_polled_I_replicate_imputed_de, is.numeric)
ttaxon_numeric_de <- ttaxon_polled_I_replicate_imputed_de[, numeric_columns_de]

##upland soils within both soils
# Set the new column as row names
ttaxon_polled_I_replicate_imputed_up <- ttaxon_polled_I_replicate_imputed %>%
  filter(soil == "up")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed_de)

# Select only numeric columns for metaMDS
numeric_columns_up <- sapply(ttaxon_polled_I_replicate_imputed_up, is.numeric)
ttaxon_numeric_up2 <- ttaxon_polled_I_replicate_imputed_up[, numeric_columns_up]



# 
# 
# 
# taxon_polled_I_replicate_up <- inner_join(taxon_I_replicate_up, taxon_poll_I_replicate_up, by = "taxon")%>%
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
# ttaxon_polled_I_replicate_up <- taxon_polled_I_replicate_up%>%
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
# ttaxon_polled_I_replicate_up$sample_for_row_names <- ttaxon_polled_I_replicate_up$sample
# 
# # Set the new column as row names
# ttaxon_polled_I_replicate_imputed_up <- ttaxon_polled_I_replicate_up %>%
#   column_to_rownames(var = "sample_for_row_names")
# 
# # Check the structure to identify non-numeric columns
# str(ttaxon_polled_I_replicate_imputed_up)
# 
# # Select only numeric columns for metaMDS
# numeric_columns_up <- sapply(ttaxon_polled_I_replicate_imputed_up, is.numeric)
# ttaxon_numeric_up <- ttaxon_polled_I_replicate_imputed_up[, numeric_columns_up]




##adonis clade I#####
#with only depressionsoil 
asv.bray1_de <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate_imputed_de[,2:(ncol(ttaxon_polled_I_replicate_imputed_de)-4)], 
                                                diag = TRUE, upper = TRUE)))

dist_matrix_I_de <- as.matrix(vegdist(ttaxon_polled_I_replicate_imputed_de[,2:(ncol(ttaxon_polled_I_replicate_imputed_de)-4)], 
                                                    diag = TRUE, upper = TRUE))


heatmap(dist_matrix_I_de)

str(dist_matrix_I_de)

adonis2(asv.bray1_de ~ clade_envde$treatment, asv.bray1_de, permutations = 9000) # p = 0.05688 . treatment
adonis2(asv.bray1_de ~ clade_envde$Cu, asv.bray1_de, permutations = 9000) # p = 0.01644 * Cu
adonis2(asv.bray1_de ~ clade_envde$waterlevel, asv.bray1_de, permutations = 9000) # p = 0.3741 waterlevel
adonis2(asv.bray1_de ~ clade_envde$waterlevel*clade_envde$Cu, asv.bray1_de, permutations = 9000) # p = 0.01489 *

adonis_I_de <- adonis2(asv.bray1_de ~ clade_envde$waterlevel+clade_envde$Cu, 
                       asv.bray1_de, permutations = 9000) # p = 0.26897 water level



adonis_I_de_df <- data.frame(
  Df = adonis_I_de$Df,
  SumOfSqs = adonis_I_de$SumOfSqs,
  R2 = adonis_I_de$R2,
  F = adonis_I_de$F,
  `Pr(>F)` = adonis_I_de$`Pr(>F)`)%>% round(2)

adonis_I_de_df <- cbind(varibles, adonis_I_de_df)

write_xlsx(adonis_I_de_df,"~/Downloads/seq/stat/adonis/adonis_I_de.xlsx")





#with only upland soils
asv.bray1_up <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate_imputed_up[,2:(ncol(ttaxon_polled_I_replicate_imputed_up)-4)], 
                                                diag = TRUE, upper = TRUE)))

dist_matrix_I_up <- as.matrix(vegdist(ttaxon_polled_I_replicate_imputed_up[,2:(ncol(ttaxon_polled_I_replicate_imputed_up)-4)], 
                                                    diag = TRUE, upper = TRUE))
heatmap(dist_matrix_I_up)


adonis2(asv.bray1_up ~ clade_env_up$Cu, asv.bray1_up, permutations = 9000) # p = 0.185 Cu
adonis2(asv.bray1_up ~ clade_env_up$waterlevel, asv.bray1_up, permutations = 9000) # p = 0.03266 * Waterlevel
adonis2(asv.bray1_up ~ clade_env_up$waterlevel*clade_env_up$Cu, asv.bray1_up, permutations = 9000) # 0.02511 * Waterlevel
adonis2(asv.bray1_up ~ clade_env_up$treatment.x, asv.bray1_up, permutations = 9000) # p = 0.0711 . treatment

adonis_I_up <- adonis2(asv.bray1_up ~ clade_env_up$waterlevel+clade_env_up$Cu, asv.bray1_up, permutations = 9000) # 0.02511 * Waterlevel


adonis_I_up_df <- data.frame(
  Df = adonis_I_up$Df,
  SumOfSqs = adonis_I_up$SumOfSqs,
  R2 = adonis_I_up$R2,
  F = adonis_I_up$F,
  `Pr(>F)` = adonis_I_up$`Pr(>F)`)%>%round(3)

adonis_I_up_df <- cbind(varibles, adonis_I_up_df)

write_xlsx(adonis_I_up_df,"~/Downloads/seq/stat/adonis/adonis_I_up.xlsx")



##adonis cladeII####
dist_matrix_II <- as.matrix(vegdist(ttaxon_polled_II_replicate_imputed[,2:(ncol(ttaxon_polled_II_replicate_imputed)-4)], 
                                    diag = TRUE, upper = TRUE))

heatmap(dist_matrix_II)

asv.bray1_II <- as.data.frame(as.matrix(vegdist(ttaxon_polled_II_replicate_imputed[,2:(ncol(ttaxon_polled_II_replicate_imputed)-4)], 
                                                diag = TRUE, upper = TRUE)))


adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$Cu, asv.bray1_II, permutations = 9000) # p = 0.2663 Cu
adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$waterlevel, asv.bray1_II, permutations = 9000) # p =  0.05155 . waterlevel
adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$waterlevel*ttaxon_polled_II_replicate_imputed$Cu,
        asv.bray1_II, permutations = 9000) #p = 0.06088  . waterlevel
adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$treatment, asv.bray1_II, permutations = 9000) # P = 0.09566 .

adonis_II <- adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$waterlevel+ttaxon_polled_II_replicate_imputed$Cu,
        asv.bray1_II, permutations = 9000) #p = 0.06088  . waterlevel


print(adonis_II)
(adonis_II$aov.tab) 
as.data.frame(adonis_II$aov.tab)
adonis_II_df <- as.data.frame(adonis_II$aov.tab)
str(adonis_II)


varibles <- c("water", "Cu", "Residual", "total") 

adonis_II_df <- data.frame(
  Df = adonis_II$Df,
  SumOfSqs = adonis_II$SumOfSqs,
  R2 = adonis_II$R2,
  F = adonis_II$F,
  `Pr(>F)` = adonis_II$`Pr(>F)`)

adonis_II_df <- cbind(varibles, adonis_II_df)

write_xlsx(adonis_II_df,"~/Downloads/seq/stat/adonis/adonis_II.xlsx")


adonis_all <- cbind(adonis_I_de_df, adonis_I_up_df, adonis_II_df)

write_xlsx(adonis_all,"~/Downloads/seq/stat/adonis/adonis_all.xlsx")


##envfit with wfps#####
##de env factors 
clade_env_de <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "de")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  mutate(Cumass = `Cu_mg/kg`)%>%
  select(-`Cu_mg/kg`, -treatment.y, -replicate, 
         -ID)%>%
  slice(-17:-32)
print(ttaxon_polled_I_replicate_imputed,  n = 32)


clade_env1_de_wf  <- clade_env_de %>%
  # slice(-17:-32) %>%
  column_to_rownames(var = "sample")%>%
  select(-amoa_Achaea, -amoa_Bacteria,-nosz_cladeI,
         -`16s_bacteria`,
         -nosz_cladeII,
         -flux,-soil, -Cu, -waterlevel)



##de cladeI NMDS####

ttaxon_numeric_II


cladeI_genus_nmds_de <- metaMDS(ttaxon_numeric_de, k = 2) #distance defualt is bray

?metaMDS
print(cladeI_genus_nmds_de)
stressplot(cladeI_genus_nmds_de)


envfit_I_de_nmds_wfps <- envfit(cladeI_genus_nmds_de, clade_env1_de_wf, permu = 999, 
                           na.rm = TRUE)


plot(envfit_I_de_nmds_wfps, display = "sites")
plot(envfit_I_de_nmds_wfps)



print(envfit_I_de_arrows_nmds)
# envfit_I_de_r_nmds <- scores(envfit_I_de_nmds$vectors$r)
# envfit_I_de_r_nmds <- sqrt(envfit_I_de_r_nmds)
# envfit_I_de_arrows_nmds$NMDS3  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS4 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1

##significnat test nmds envfit cladeI de####
r_I_de_nmds_env_wfps <- as.matrix(envfit_I_de_nmds_wfps$vectors$r) 
p_I_de_nmds_env_wfps <- as.matrix(envfit_I_de_nmds_wfps$vectors$pvals)
env_p_I_de_nmds_wfps <- cbind(r_I_de_nmds_env_wfps, p_I_de_nmds_env_wfps)
colnames(env_p_I_de_nmds_wfps) <- c("r2","p_value")
K_I_de_nmds_wfps <- as.data.frame(env_p_I_de_nmds_wfps)%>%
  mutate(p_adjust = p.adjust(K_I_de_nmds_wfps$p_value, method = "BH"))
K_I_de_nmds_wfps

write_xlsx(K_I_de_nmds_wfps,"~/Downloads/seq/stat/envfit/wfps/K_I_de_nmds_wfps2.xlsx")




envfit_I_de_arrows_nmds_wfps <- as.data.frame(envfit_I_de_nmds_wfps$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_de_arrows_nmds_wfps <- envfit_I_de_arrows_nmds_wfps%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



##envfit cladeI_genus_nmds_de, and plotting######
#sp
cladeI_genus_nmds_points_de_sp <- as.data.frame(scores(cladeI_genus_nmds_de)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_sp1 <- cladeI_genus_nmds_points_de_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_de_si <- as.data.frame(scores(cladeI_genus_nmds_de)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_si1 <- cladeI_genus_nmds_points_de_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  

centroid_I_de2 <- cladeI_genus_nmds_points_de%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )


######de NMDs species plotting +envfit #####
cladeI_genus_nmds_points_de_sp%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2)) +
  # geom_point(data = centroid_I_de2, size = 5, show.legend = FALSE)+
  geom_segment(data = cladeI_genus_nmds_points_de_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_de_arrows_nmds_wfps,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_de_arrows_nmds_wfps,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_de_sp1,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE,
                   size = 3.5
  ) +
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de_wfps, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/2/cladeI_de_nmds_env_species_wfps.tiff",  width=7, height=8)
ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/2/cladeI_de_nmds_env_species_wpfs.jpeg",  width=7, height=8)




##upland environmental varibels######
clade_env_up <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "up")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  mutate(Cumass = `Cu_mg/kg`)%>%
  select(-`Cu_mg/kg`, -treatment.y, -replicate, 
         -ID,  -amoa_Achaea, -nosz_cladeII)%>%
  slice(-14:-32) 


clade_env1_up_wf <- clade_env_up %>%
  # slice(-14:-32) %>%
  column_to_rownames(var = "sample")%>%
  select( -amoa_Bacteria,-nosz_cladeI,
          -`16s_bacteria`,
          -flux,-soil, -Cu, -waterlevel,)



##up cladeI NMDS####
cladeI_genus_nmds_up <- metaMDS(ttaxon_numeric_up, k = 2) #distance defualt is bray

?metaMDS
plot(cladeI_genus_nmds_up)



envfit_I_up_nmds_wfps <- envfit(cladeI_genus_nmds_up, clade_env1_up_wf, permu = 999, 
                                na.rm = TRUE)


plot(envfit_I_up_nmds_wfps, display = "sites")
plot(envfit_I_up_nmds_wfps)



envfit_I_up_arrows_nmds_wfps <- as.data.frame(envfit_I_up_nmds_wfps$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_up_arrows_nmds_wfps <- envfit_I_up_arrows_nmds_wfps%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_I_up_arrows_nmds)
# envfit_I_up_r_nmds <- scores(envfit_I_up_nmds$vectors$r)
# envfit_I_up_r_nmds <- sqrt(envfit_I_up_r_nmds)
# envfit_I_up_arrows_nmds$NMDS3  <- envfit_I_up_arrows_nmds$NMDS1  / as.numeric(envfit_I_up_r_nmds) * 1
# envfit_I_up_arrows_nmds$NMDS4 <- envfit_I_up_arrows_nmds$NMDS2 / as.numeric(envfit_I_up_r_nmds) * 1

##significnat test nmds envfit cladeI up####
r_I_up_nmds_env_wfps <- as.matrix(envfit_I_up_nmds_wfps$vectors$r) 
p_I_up_nmds_env_wfps <- as.matrix(envfit_I_up_nmds_wfps$vectors$pvals)
env_p_I_up_nmds_wfps <- cbind(r_I_up_nmds_env_wfps, p_I_up_nmds_env_wfps)
colnames(env_p_I_up_nmds_wfps) <- c("r2","p_value")
K_I_up_nmds_wfps <- as.data.frame(env_p_I_up_nmds_wfps)%>%
  mutate(p_adjust = p.adjust(K_I_up_nmds_wfps$p_value, method = "BH"))
K_I_up_nmds_wfps

write_xlsx(K_I_up_nmds_wfps,"~/Downloads/seq/stat/envfit/wfps/K_I_up_nmds_wfps2.xlsx")


##envfit cladeI_genus_nmds_up, and plotting######
#sp
cladeI_genus_nmds_points_up_sp <- as.data.frame(scores(cladeI_genus_nmds_up)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_sp1 <- cladeI_genus_nmds_points_up_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_up_si <- as.data.frame(scores(cladeI_genus_nmds_up)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_si1 <- cladeI_genus_nmds_points_up_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  

centroid_I_up2 <- cladeI_genus_nmds_points_up%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )


######up NMDs species plotting +envfit #####
cladeI_genus_nmds_points_up_sp%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2)) +
  geom_segment(data = cladeI_genus_nmds_points_up_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_up_arrows_nmds_wfps,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_up_arrows_nmds_wfps,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_up_sp1,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE,
                   size = 3.5
  ) +
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Upland soils, using upland soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_up, clade_env1_up_wfps, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/2/cladeI_up_nmds_env_species_wfps.tiff",  width=7, height=8)
ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/2/cladeI_up_nmds_env_species_wpfs.jpeg",  width=7, height=8)



####clade II#####

cladeII_genus_nmds <- metaMDS(ttaxon_numeric_II, k = 2)

envfit_II_de_nmds_wfps <- envfit(cladeII_genus_nmds, clade_env1_de_wf, permu = 999, 
                            na.rm = TRUE)

r_II_de_nmds_env_wfps <- as.matrix(envfit_II_de_nmds_wfps$vectors$r) 
p_II_de_nmds_env_wfps <- as.matrix(envfit_II_de_nmds_wfps$vectors$pvals)
env_p_II_de_nmds_wfps <- cbind(r_II_de_nmds_env_wfps,p_II_de_nmds_env_wfps)
colnames(env_p_II_de_nmds_wfps) <- c("r2","p_value")
K_II_de_nmds_wfps <- as.data.frame(env_p_II_de_nmds_wfps)%>%
  mutate(p_adjust = p.adjust(K_II_de_nmds_wfps$p_value, method = "BH"))
K_II_de_nmds_wfps

write_xlsx(K_II_de_nmds_wfps,"~/Downloads/seq/stat/envfit/wfps/K_II_de_nmds_wfps2.xlsx")




envfit_II_de_arrows_nmds_wfps <- as.data.frame(envfit_II_de_nmds_wfps$vectors$arrows) %>%
  mutate(varNames = rownames(.))
envfit_II_de_arrows_nmds_wfps <- envfit_II_de_arrows_nmds_wfps%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_II_de_arrows_nmds_wfps)
envfit_II_de_r_nmds_wfps <- scores(envfit_II_de_nmds_wfps$vectors$r)
envfit_II_de_r_nmds_wfps <- sqrt(envfit_II_de_r_nmds_wfps)
# envfit_I_de_arrows_nmds$NMDS1  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS2 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1





##envfit cladeII_genus_nmds_de, and plotting######
#sp
cladeII_genus_nmds_points_de_sp <- as.data.frame(scores(cladeII_genus_nmds)$species)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_sp1 <- cladeII_genus_nmds_points_de_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_(.*)", "\\1 \\2"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeII_genus_nmds_points_de_si <- as.data.frame(scores(cladeII_genus_nmds)$site)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_si1 <- cladeII_genus_nmds_points_de_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  



######cladeII de NMDs species plotting +envfit #####
cladeII_genus_nmds_points_de_sp%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2)) +
  geom_segment(data = cladeII_genus_nmds_points_de_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_II_de_arrows_nmds_wfps,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_II_de_arrows_nmds_wfps,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeII_genus_nmds_points_de_sp1,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE,
                   size = 3.5
  ) +
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(title = "*nosZ* clade II",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeII_genus_nmds_de, clade_env1_de_wfps, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/cladeII_de_nmds_env_species_wfps.tiff",  width=7, height=8)
ggsave("~/Downloads/seq/envfitplot_dot/Wwfps/cladeII_de_nmds_env_species_wpfs.jpeg",  width=7, height=8)



#with flux -wfps#####
##de env factors 
clade_env_de <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "de")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  mutate(Cumass = `Cu_mg/kg`)%>%
  select(-`Cu_mg/kg`, -treatment.y, -replicate, -treatment.x,
         -ID, -soil, -Cu, -waterlevel)


clade_env1_de_f  <- clade_env_de %>%
  slice(-17:-32) %>%
  column_to_rownames(var = "sample")%>%
  select(-amoa_Achaea, -amoa_Bacteria,-nosz_cladeI,
         -`16s_bacteria`,
         -nosz_cladeII,
         -wfps)



##de cladeI NMDS####
cladeI_genus_nmds_de <- metaMDS(ttaxon_numeric_de, k = 2) #distance defualt is bray

?metaMDS
plot(cladeI_genus_nmds_de)



envfit_I_de_nmds_f <- envfit(cladeI_genus_nmds_de, clade_env1_de_f, permu = 999, 
                                na.rm = TRUE)


plot(envfit_I_de_nmds_f, display = "sites")
plot(envfit_I_de_nmds_f)



envfit_I_de_arrows_nmds_f <- as.data.frame(envfit_I_de_nmds_f$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_de_arrows_nmds_f <- envfit_I_de_arrows_nmds_f%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))




print(envfit_I_de_arrows_nmds)
# envfit_I_de_r_nmds <- scores(envfit_I_de_nmds$vectors$r)
# envfit_I_de_r_nmds <- sqrt(envfit_I_de_r_nmds)
# envfit_I_de_arrows_nmds$NMDS3  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS4 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1

##significnat test nmds envfit cladeI de####
r_I_de_nmds_env_f <- as.matrix(envfit_I_de_nmds_f$vectors$r) 
p_I_de_nmds_env_f <- as.matrix(envfit_I_de_nmds_f$vectors$pvals)
env_p_I_de_nmds_f <- cbind(r_I_de_nmds_env_f, p_I_de_nmds_env_f)
colnames(env_p_I_de_nmds_f) <- c("r2","p_value")
K_I_de_nmds_f <- as.data.frame(env_p_I_de_nmds_f)%>%
  mutate(p_adjust = p.adjust(K_I_de_nmds_f$p_value, method = "BH"))
K_I_de_nmds_f

write_xlsx(K_I_de_nmds_f,"~/Downloads/seq/stat/envfit/flux/K_I_de_nmds_f2.xlsx")




##envfit cladeI_genus_nmds_de, and plotting######
#sp
cladeI_genus_nmds_points_de_sp <- as.data.frame(scores(cladeI_genus_nmds_de)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_sp1 <- cladeI_genus_nmds_points_de_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_de_si <- as.data.frame(scores(cladeI_genus_nmds_de)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_de_si1 <- cladeI_genus_nmds_points_de_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  

centroid_I_de2 <- cladeI_genus_nmds_points_de%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )





######de NMDs flux species plotting +envfit #####
cladeI_genus_nmds_points_de%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  geom_point(data = centroid_I_de2, size = 5, show.legend = FALSE)+
  geom_segment(data = cladeI_genus_nmds_points_de_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_de_arrows_nmds_f,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_de_arrows_nmds_f,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_de_sp1,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.8,parse = TRUE,
                   size = 3,
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
                     labels=c("Depression soil Cu0 WHC60",
                              "Depression soil Cu260 WHC60",
                              "Depression soil Cu0 WHC90",
                              "Depression soil Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de_f, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeI_de_nmds_env_species_f2.tiff",  width=9, height=9)
ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeI_de_nmds_env_species_f2.jpeg",  width=9, height=9)




##upland environmental varibels
clade_env_up <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")%>%
  filter(soil == "up")%>%
  mutate(ID = str_replace(sample, "\\w(\\d*)", "\\1"))%>%
  select(-soil, -Cu, -waterlevel)%>%
  right_join(env, by = "ID")%>%
  mutate(Cumass = `Cu_mg/kg`)%>%
  select(-`Cu_mg/kg`, -treatment.y, -replicate, -treatment.x,
         -ID, -soil, -Cu, -waterlevel, -amoa_Achaea, -nosz_cladeII)


clade_env1_up_f <- clade_env_up %>%
  slice(-14:-32) %>%
  column_to_rownames(var = "sample")%>%
  select( -amoa_Bacteria,-nosz_cladeI,
          -`16s_bacteria`,
          -wfps)



##up cladeI NMDS####
cladeI_genus_nmds_up <- metaMDS(ttaxon_numeric_up, k = 2) #distance defualt is bray

?metaMDS
plot(cladeI_genus_nmds_up)



envfit_I_up_nmds_f <- envfit(cladeI_genus_nmds_up, clade_env1_up_f, permu = 999, 
                                na.rm = TRUE)


plot(envfit_I_up_nmds_f, display = "sites")
plot(envfit_I_up_nmds_f)





envfit_I_up_arrows_nmds_f <- as.data.frame(envfit_I_up_nmds_f$vectors$arrows) %>%
  mutate(varNames = rownames(.))
envfit_I_up_arrows_nmds_f <- envfit_I_up_arrows_nmds_f%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_I_up_arrows_nmds)
# envfit_I_up_r_nmds <- scores(envfit_I_up_nmds$vectors$r)
# envfit_I_up_r_nmds <- sqrt(envfit_I_up_r_nmds)
# envfit_I_up_arrows_nmds$NMDS3  <- envfit_I_up_arrows_nmds$NMDS1  / as.numeric(envfit_I_up_r_nmds) * 1
# envfit_I_up_arrows_nmds$NMDS4 <- envfit_I_up_arrows_nmds$NMDS2 / as.numeric(envfit_I_up_r_nmds) * 1

##significnat test nmds envfit cladeI de####
r_I_up_nmds_env_f <- as.matrix(envfit_I_up_nmds_f$vectors$r) 
p_I_up_nmds_env_f <- as.matrix(envfit_I_up_nmds_f$vectors$pvals)
env_p_I_up_nmds_f <- cbind(r_I_up_nmds_env_f, p_I_up_nmds_env_f)
colnames(env_p_I_up_nmds_f) <- c("r2","p_value")
K_I_up_nmds_f <- as.data.frame(env_p_I_up_nmds_f)%>%
  mutate(p_adjust = p.adjust(K_I_up_nmds_f$p_value, method = "BH"))
K_I_up_nmds_f

write_xlsx(K_I_up_nmds_f,"~/Downloads/seq/stat/envfit/flux/K_I_up_nmds_f2.xlsx")



##envfit cladeI_genus_nmds_up, and plotting######
#sp
cladeI_genus_nmds_points_up_sp <- as.data.frame(scores(cladeI_genus_nmds_up)$species)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_sp1 <- cladeI_genus_nmds_points_up_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeI_genus_nmds_points_up_si <- as.data.frame(scores(cladeI_genus_nmds_up)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_si1 <- cladeI_genus_nmds_points_up_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  

centroid_I_up2 <- cladeI_genus_nmds_points_up%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )


######up NMDs species plotting +envfit #####
# cladeI_genus_nmds_points_up_sp%>%
  cladeI_genus_nmds_points_up%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  geom_point(data = centroid_I_up, size = 5, show.legend = FALSE)+
  geom_segment(data = cladeI_genus_nmds_points_up_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_I_up_arrows_nmds_f,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_I_up_arrows_nmds_f,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_up_sp1,
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
  scale_color_manual(name = NULL, values  = c(upland0Cu60WHC_c1,
                                              upland260Cu60WHC_c1,
                                              upland0Cu90WHC_c1,
                                              upland260Cu90WHC_c1),
                     labels=c("upland soil Cu0 WHC60",
                              "upland soil Cu260 WHC60",
                              "upland soil Cu0 WHC90",
                              "upland soil Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(upland0Cu60WHC_c1,
                                             upland260Cu60WHC_c1,
                                             upland0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Upland soils, using upland soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_up, clade_env1_up_f, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeI_up_nmds_env_species_f2.tiff",  width=9, height=9)
ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeI_up_nmds_env_species_f2.jpeg",  width=9, height=9)



####with flux clade II#####


cladeII_genus_nmds <- metaMDS(ttaxon_numeric_II, k = 2)

envfit_II_de_nmds_f <- envfit(cladeII_genus_nmds, clade_env1_de_f, permu = 999, 
                                 na.rm = TRUE)

r_II_de_nmds_env_f <- as.matrix(envfit_II_de_nmds_f$vectors$r) 
p_II_de_nmds_env_f <- as.matrix(envfit_II_de_nmds_f$vectors$pvals)
env_p_II_de_nmds_f <- cbind(r_II_de_nmds_env_f,p_II_de_nmds_env_f)
colnames(env_p_II_de_nmds_f) <- c("r2","p_value")
K_II_de_nmds_f <- as.data.frame(env_p_II_de_nmds_f)%>%
  mutate(p_adjust = p.adjust(K_II_de_nmds_f$p_value, method = "BH"))
K_II_de_nmds_f

write_xlsx(K_II_de_nmds_f,"~/Downloads/seq/stat/envfit/flux/K_II_de_nmds_f2.xlsx")




envfit_II_de_arrows_nmds_f <- as.data.frame(envfit_II_de_nmds_f$vectors$arrows) %>%
  mutate(varNames = rownames(.))
envfit_II_de_arrows_nmds_f <- envfit_II_de_arrows_nmds_f%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



print(envfit_II_de_arrows_nmds_f)
envfit_II_de_r_nmds_f <- scores(envfit_II_de_nmds_f$vectors$r)
envfit_II_de_r_nmds_f <- sqrt(envfit_II_de_r_nmds_f)
# envfit_I_de_arrows_nmds$NMDS1  <- envfit_I_de_arrows_nmds$NMDS1  / as.numeric(envfit_I_de_r_nmds) * 1
# envfit_I_de_arrows_nmds$NMDS2 <- envfit_I_de_arrows_nmds$NMDS2 / as.numeric(envfit_I_de_r_nmds) * 1





##envfit cladeII_genus_nmds_de, and plotting######
#sp
cladeII_genus_nmds_points_de_sp <- as.data.frame(scores(cladeII_genus_nmds)$species)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_sp1 <- cladeII_genus_nmds_points_de_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_(.*)", "\\1 \\2"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeII_genus_nmds_points_de_si <- as.data.frame(scores(cladeII_genus_nmds)$site)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_si1 <- cladeII_genus_nmds_points_de_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  



######de NMDs species plotting +envfit #####
# cladeII_genus_nmds_points_de_sp%>%
  cladeII_genus_nmds_points%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color = treatment)) +
  geom_point(data = centroid_II_de, size = 5, show.legend = FALSE)+
  geom_segment(data = cladeII_genus_nmds_points_de_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7
  ) +
  geom_segment(data = envfit_II_de_arrows_nmds_f,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color="#B2182B") +
  geom_label_repel(data = envfit_II_de_arrows_nmds_f,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#B2182B",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeII_genus_nmds_points_de_sp1,
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
                     labels=c("Depression soil Cu0 WHC60",
                              "Depression soil Cu260 WHC60",
                              "Depression soil Cu0 WHC90",
                              "Depression soil Cu260 WHC90"))+
  scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  # geom_vline(aes(xintercept = 0), linetype = "dotted") +
  # geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(title = "*nosZ* clade II",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeII_genus_nmds_de, clade_env1_de_f, permu = 999, na.rm = TRUE)") +
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


ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeII_de_nmds_env_species_flux2.tiff",  width=7, height=8)
ggsave("~/Downloads/seq/envfitplot_dot/Wflux/cladeII_de_nmds_env_species_flux2.jpeg",  width=7, height=8)





##simper####
#de cladeI 
sim_I_de <- simper(ttaxon_numeric_de, clade_env_de$treatment.x)




summary(sim_I_de)$dell_dehl %>%
  # round(3) %>%
  head()

comparisons_de <- c("dell_dehl", "dell_delh", "dell_dehh", "dehl_delh","dehl_dehh", "delh_dehh")

simper_results_I_de <- c()

for(i in 1:length(comparisons_de)) {
  require(tidyverse)
  temp <- summary(sim_I_de)[as.character(comparisons_de[i])] %>%
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons_de[i],".", sep = ""), "", colnames(temp))
  temp <- temp %>%
    mutate(Comparison = comparisons_de[i],
           Position = row_number()) %>%
    rownames_to_column(var = "Species")
  simper_results_I_de <- rbind(simper_results_I_de, temp)
}



simper_results_I_de


simper_results_I_de%>%group_by(Comparison) %>%
  summarize(sum.average = sum(average))


simper_results_I_de%>%
  filter(Position == "1")


ggplot(data = simper_results_I_de,
             aes(x = Position, y = cumsum)) +
  geom_line(aes(colour = Comparison)) +
  theme_bw()

ggplot(data = simper_results_I_de,
             aes(x = Position, y = average)) +
  geom_line(aes(colour = Comparison)) +
  theme_bw()


###de by water level

sim_I_de_wl <- simper(ttaxon_numeric_de2, clade_env_de$waterlevel)

summary(sim_I_de_wl)
sim_I_de_wl$"60_90"





###de by Cu

sim_I_de_Cu <- simper(ttaxon_numeric_de2, clade_env_de$Cu)

summary(sim_I_de_Cu)

sim_I_de_Cu$"0_260"




#up cladeI 
sim_I_up <- simper(ttaxon_numeric_up, clade_env_up$treatment.x)


summary(sim_I_up)$upll_uphl %>%
  # round(3) %>%
  head()

comparisons_up <- c("upll_uphl", "upll_uplh", "upll_uphh", "uphl_uplh","uphl_uphh", "uplh_uphh")


simper_results_I_up <- c()

for(i in 1:length(comparisons_up)) {
  require(tidyverse)
  temp <- summary(sim_I_up)[as.character(comparisons_up[i])] %>%
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons_up[i],".", sep = ""), "", colnames(temp))
  temp <- temp %>%
    mutate(Comparison = comparisons_up[i],
           Position = row_number()) %>%
    rownames_to_column(var = "Species")
  simper_results_I_up <- rbind(simper_results_I_up, temp)
}



simper_results_I_up


simper_results_I_up%>%group_by(Comparison) %>%
  summarize(sum.average = sum(average))

simper_results_I_up%>%
  filter(Position == "1")




###up by water level

sim_I_up_wl <- simper(ttaxon_numeric_up, clade_env_up$waterlevel)

summary(sim_I_up_wl)
sim_I_up_wl$"60_90"





###up by Cu

sim_I_up_Cu <- simper(ttaxon_numeric_up, clade_env_up$Cu)

summary(sim_I_up_Cu)
summary(sim_I_de_Cu)
sim_I_up_Cu$"0_260"
summary(sim_I_de_wl)
summary(sim_I_up_wl)







###clade II
sim_II_de <- simper(ttaxon_numeric_II, clade_env_de$treatment.x)


summary(sim_II_de)$dell_dehl %>%
  # round(3) %>%
  head()

comparisons_de <- c("dell_dehl", "dell_delh", "dell_dehh", "dehl_delh","dehl_dehh", "delh_dehh")

simper_results_II_de <- c()

for(i in 1:length(comparisons_de)) {
  require(tidyverse)
  temp <- summary(sim_II_de)[as.character(comparisons_de[i])] %>%
    as.data.frame()
  colnames(temp) <- gsub(
    paste(comparisons_de[i],".", sep = ""), "", colnames(temp))
  temp <- temp %>%
    mutate(Comparison = comparisons_de[i],
           Position = row_number()) %>%
    rownames_to_column(var = "Species")
  simper_results_II_de <- rbind(simper_results_II_de, temp)
}



simper_results_II_de


simper_results_II_de%>%group_by(Comparison) %>%
  summarize(sum.average = sum(average))


simper_results_II_de%>%
  filter(Position == "1")
