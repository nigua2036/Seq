source("script/packages.R")

#read in data file#####
cI <-  read.table("~/Downloads/newn2o/total-spruce-taxa_clade_I_filtered.txt",
                  header=TRUE, sep="\t")

cI <- cI  %>%
  dplyr::select(taxa,clade,value,sample, soil, Cu, waterlevel)



# cI_filtered <-  cI%>%
#   subset(sample == "n49")
# 
# write.table(cI_filtered, "~/Downloads/newn2o/total-spruce-taxa_clade_I_filtered.txt")

cI <- cI%>%
  mutate(treatment =  paste0(cI$soil,cI$Cu, cI$waterlevel))%>%
  group_by(sample)%>%
  mutate(rel_abund_raw = (value / sum(value))*100)  


# Cu_gene_f <- Cu_gene%>%
#   subset(ID == "49")
# 
# write_xlsx(Cu_gene_f, "~/Downloads/newn2o/Cu_gene_filtered.xlsx")

##copy number from qPCR
Cu_gene <- read_excel("~/Downloads/newn2o/Cu_gene_filtered.xlsx")

cpclade <- Cu_gene%>%
  mutate(sample = as.character(ID),
         cladeI = nosz_cladeI,
         cladeII = nosz_cladeII )%>%
  mutate(sample = str_c("n", sample))%>%
  dplyr::select(sample, cladeI, cladeII)%>%
  mutate(cladeI = as.numeric(cladeI),
         cladeII = as.numeric(cladeII))


cI <- cI%>%
  left_join(cpclade, by = "sample")

cI <- cI%>%
  mutate(cladeIcp_raw = rel_abund_raw*cladeI ) #cladeIcp_raw copy number for each taxa





meancI <-  cI %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100,
         N = N)%>% #rel_abund relative abundnce including NAs
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund,N, treatment)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
               names_to="level",
               values_to="taxon")


# meancI_wide same as meancI but different format
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



#taxa,genus plotting#####

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
  # mutate(replicates = case_when(
  #   treatment == "dehh" ~ 4,
  #   treatment == "delh" ~ 4,
  #   treatment == "dell" ~ 4,
  #   treatment == "dehl" ~ 4,
  #   treatment == "uplh" ~ 3,
  #   treatment == "upll" ~ 4,
  #   treatment == "uphl" ~ 2,
  #   treatment == "uphh" ~ 4
  #   ))%>% #when you have different replicate numbers, here only sample 49, one replicate
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))



taxon_pollIg2 <- taxon_rel_abundI_genus2%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3, #pool relative abundance < 3% into Other
            mean = mean(mean_rel_abund3),
            .groups = "drop")

#set color
library(RColorBrewer)

colors_set3 <- brewer.pal(12, "Set3")
colors_paired <- brewer.pal(1, "Paired")  # Adding one color from another palette

# Combine the colors to have a total of 13
all_colors <- c(colors_set3, colors_paired)
all_colors <- c("#8DD3C7", "#FFFFB3",  "#1F78B4",  "#A6CEE3", "#FB8072", "#7570B3",
                "#FDB462" ,"#B3DE69" ,"#FCCDE5", "#BC80BD" , "#FFED6F", "#E5C494","gray","#FED9A6")


#change treatment into your treatment
inner_join(taxon_rel_abundI_genus2, taxon_pollIg2, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  # mutate(treatment = factor(treatment,
  #                           levels = 
  #                             c("upll",
  #                               "uphl",
  #                               "uplh",
  #                               "uphh",
  #                               "dell",
  #                               "dehl",
  #                               "delh",
  #                               "dehh")))%>%
  group_by(treatment, taxon)%>%
  summarise(mean_rel_abund = sum(mean_rel_abund3),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon,mean, .desc = TRUE), 
         taxon = fct_shift(taxon, n = 1)) %>%
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
  # scale_x_discrete(breaks=c("dehh",
  #                           "dehl",
  #                           "delh",
  #                           "dell",
  #                           "uphh",
  #                           "uphl",
  #                           "uplh",
  #                           "upll"),
  #                  labels = c("depression<br>Cu260<br>water90",
  #                             "depression<br>Cu260<br>water60",
  #                             "depression<br>Cu0<br>water90",
  #                             "depression<br>Cu0<br>water60",
  #                             "upland<br>Cu260<br>water90",
  #                             "upland<br>Cu260<br>water60",
  #                             "upland<br>Cu0<br>water90",
  #                             "upland<br>Cu0<br>water60"))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        axis.text.y = element_markdown(size = 12),
        legend.text = element_markdown(size = 12),
        legend.key.size = unit(28, "pt"),
        axis.title.y  = element_markdown(size = 12),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))





#NMDs species plotting  + envfit#####
# from Riffomonas video
# choose your level, e.g., genus here for later analysis

taxon_rel_abundI_genus3 <- meancI %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>% #rel_abund2 relarive abundance without NA
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
  # mutate(replicates = case_when(
  #   treatment == "dehh" ~ 4,
  #   treatment == "delh" ~ 4,
  #   treatment == "dell" ~ 4,
  #   treatment == "dehl" ~ 4,
  #   treatment == "uplh" ~ 3,
  #   treatment == "upll" ~ 4,
  #   treatment == "uphl" ~ 2,
  #   treatment == "uphh" ~ 4))%>%
  mutate(mean_rel_abund3 = mean_rel_abund2/replicates)




taxon_pollIg3 <- taxon_rel_abundI_genus3%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3, #pool relative abundance <3% into Other
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





# my treatment list
tlist <- meancI_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)




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






# environmental data for later
microbe <- read_excel("~/Downloads/newn2o/microbe_filtered.xlsx")

# microbe_filtered <- microbe%>%
#   subset(ID == "49")
#   
# microbe_filtered 
# write_xlsx(microbe_filtered, "~/Downloads/newn2o/microbe_filtered.xlsx")

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






ttaxon_polled_I_replicate2_imputed$sample_for_row_names <- ttaxon_polled_I_replicate2_imputed$sample

# Set the new column as row names
ttaxon_polled_I_replicate2_imputed <- ttaxon_polled_I_replicate2_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate2_imputed)

# Select only numeric columns for metaMDS
numeric_columns <- sapply(ttaxon_polled_I_replicate2_imputed, is.numeric)
ttaxon_numeric2 <- ttaxon_polled_I_replicate2_imputed[, numeric_columns]


# Set the new column as row names
ttaxon_polled_I_replicate2_imputed_de <- ttaxon_polled_I_replicate2_imputed %>%
  filter(soil == "de")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate2_imputed_de)

# Select only numeric columns for metaMDS
numeric_columns_de <- sapply(ttaxon_polled_I_replicate2_imputed_de, is.numeric)
ttaxon_numeric2_de2 <- ttaxon_polled_I_replicate2_imputed_de[, numeric_columns_de]





##adonis clade I#####
#with only depressionsoil 
asv.bray1_de2 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_de[,2:(ncol(ttaxon_polled_I_replicate2_imputed_de)-4)], 
                                                 diag = TRUE, upper = TRUE)))

dist_matrix_I_de2 <- as.matrix(vegdist(ttaxon_polled_I_replicate2_imputed_de[,2:(ncol(ttaxon_polled_I_replicate2_imputed_de)-4)], 
                                       diag = TRUE, upper = TRUE))


heatmap(dist_matrix_I_de2)

str(dist_matrix_I_de2)

adonis2(asv.bray1_de2 ~ clade_env_de$treatment.x, asv.bray1_de2, permutations = 9000) #  treatment
adonis2(asv.bray1_de2 ~ clade_env_de$Cu, asv.bray1_de2, permutations = 9000) #  Cu
adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel, asv.bray1_de2, permutations = 9000) #  waterlevel
adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel*clade_env_de$Cu, asv.bray1_de2, permutations = 9000) 

adonis_I_de2 <- adonis2(asv.bray1_de2 ~ clade_env_de$waterlevel+clade_env_de$Cu, 
                        asv.bray1_de2, permutations = 9000) #  water level



adonis_I_de_df2 <- data.frame(
  Df = adonis_I_de2$Df,
  SumOfSqs = adonis_I_de2$SumOfSqs,
  R2 = adonis_I_de2$R2,
  F = adonis_I_de2$F,
  `Pr(>F)` = adonis_I_de2$`Pr(>F)`)%>% round(2)

varibles <- c("water", "Cu", "Residual", "total") 

adonis_I_de_df2 <- cbind(varibles, adonis_I_de_df2)

 write_xlsx(adonis_I_de_df2,"~/Downloads/seq/stat/adonis2/adonis_I_de2.xlsx")





##*de cladeI NMDS####
cladeI_genus_nmds_de2 <- metaMDS(ttaxon_numeric2_de2, k = 2) #distance default is bray

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
  # mutate(treatment = factor(treatment,
  #                           levels = 
  #                             c(
  #                               "dell",
  #                               "dehl",
  #                               "delh",
  #                               "dehh")))%>%
  # slice(-17:-29)

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
  
  
  # scale_color_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
  #                                             depression260Cu60WHC_c1,
  #                                             depression0Cu90WHC_c1,
  #                                             depression260Cu90WHC_c1),
  #                    labels=c("Cu0 WHC60",
  #                             "Cu260 WHC60",
  #                             "Cu0 WHC90",
  #                             "Cu260 WHC90"))+
  # scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
  #                                            depression260Cu60WHC_c1,
  #                                            depression0Cu90WHC_c1,
  #                                            depression260Cu90WHC_c1),
  #                   guide = "none")+
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





