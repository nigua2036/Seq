
###cladeI #####
#read in
cI <-  read.table("~/Downloads/seq/total-spruce-taxa_clade_I.txt",
                  header=TRUE, sep="\t")%>%
  dplyr::select(taxa,clade,value,sample, soil, Cu, waterlevel)%>%
  mutate(treatment =  paste0(cI$soil,cI$Cu, cI$waterlevel))%>%
  group_by(sample)%>%
  mutate(rel_abund = (value / sum(value))*100)


#not real mean, bad naming
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


###taxon with replicates
taxon_I_replicate  <- meancI %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 

#ready to poll
taxon_poll_I_replicate <- taxon_I_replicate%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")


#poll rel_abund2 <3 (%) into "Other"
taxon_polled_I_replicate <- inner_join(taxon_I_replicate, taxon_poll_I_replicate, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(sample, taxon)%>%
  summarise(mean_rel_abund = sum(rel_abund2),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1))

#adding sample list
tlist <- meancI_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)


ttaxon_polled_I_replicate <- taxon_polled_I_replicate%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist, by = "sample") 



# Permanova test: are samples significantly different from each other 

library("vegan")

#replace NAs with 0
ttaxon_polled_I_replicate_imputed <- ttaxon_polled_I_replicate %>%
  replace_na(list(" Mesorhizobium" = 0, " Microvirga"  = 0,
                  " Noviherbaspirillum" = 0, " Polaromonas" = 0,
                  " Pseudomonas"   = 0,    " Ramlibacter"  = 0,     " Rhizobium/Agrobacterium_group" = 0, " Skermanella"    = 0,
                  "Other"  = 0,  " unclassified_Acetobacteraceae"  = 0, " Azospirillum"   = 0,  " Neisseria" = 0 ))




###adonis test
asv.bray1 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate_imputed[,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)], 
                                             diag = TRUE, upper = TRUE)))

adonis2(asv.bray1 ~ ttaxon_polled_I_replicate_imputed$soil, asv.bray1, permutations = 3000) # p = 0.0003332 ***
adonis2(asv.bray1 ~ ttaxon_polled_I_replicate_imputed$Cu, asv.bray1, permutations = 3000) # p = 0.3176
adonis2(asv.bray1 ~ ttaxon_polled_I_replicate_imputed$waterlevel, asv.bray1, permutations = 3000) # p = 0.05631 .
adonis2(asv.bray1 ~ ttaxon_polled_I_replicate_imputed$treatment, asv.bray1, permutations = 3000) # p = 0.0009997 ***

adonis2(asv.bray1 ~ ttaxon_polled_I_replicate_imputed$soil, asv.bray1, permutations = 3000, 
        strata = ttaxon_polled_I_replicate_imputed$waterlevel) # p = 0.0003332 ***


## Run the NMDS####
###cladeI_genus_nmds with relative abundance####
# Create a new column for row names
ttaxon_polled_I_replicate_imputed$sample_for_row_names <- ttaxon_polled_I_replicate_imputed$sample

# Set the new column as row names
ttaxon_polled_I_replicate_imputed <- ttaxon_polled_I_replicate_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed)

# Select only numeric columns for metaMDS
numeric_columns <- sapply(ttaxon_polled_I_replicate_imputed, is.numeric)
ttaxon_numeric <- ttaxon_polled_I_replicate_imputed[, numeric_columns]

# Run metaMDS on the numeric columns
cladeI_genus_nmds <- metaMDS(ttaxon_numeric, k = 2)

# View the results
print(cladeI_genus_nmds)
stressplot(cladeI_genus_nmds)
#stress = 0.09530978 


###another way for metaMDS cladeI_genus_nmds2####
# asv.bray1 <- as.data.frame(as.matrix(vegdist(ttaxon_polled_I_replicate_imputed
#                                              [,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)], 
#                                              diag = TRUE, upper = TRUE)))
cladeI_genus_nmds2 <- metaMDS(asv.bray1,  k = 2)
stressplot(cladeI_genus_nmds2)
print(cladeI_genus_nmds2)
# stress = 0.0884576 


#cladeI_genus_nmds from metaMDS(relative abundance)
cladeI_genus_nmds_points <- as.data.frame(cladeI_genus_nmds$points)%>%
  mutate(sample = row.names(.))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c("upll",
                                "uphl",
                                "uplh",
                                "uphh",
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))


#cladeI_genus_nmds2 from metaMDS(vegdist(relative abundance))
#I wonder what is the difference between these 2?
cladeI_genus_nmds_points2 <- as.data.frame(cladeI_genus_nmds2$points)%>%
  mutate(sample = row.names(cladeI_genus_nmds$points))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c("upll",
                                "uphl",
                                "uplh",
                                "uphh",
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))




# Stress plot
stressplot(cladeI_genus_nmds)

# svg("plots/asv_nmds_stressplot.svg", width = 4, height = 4)
stressplot(cladeI_genus_nmds)
dev.off()

stressplot(cladeI_genus_nmds2)


##define colors####
upland0Cu60WHC_c1 <- "#BC80BD"

upland0Cu90WHC_c1 <-  "#E6AB02"

upland260Cu60WHC_c1 <-  "#377EB8"

upland260Cu90WHC_c1  <-   "#66A61E"


depression0Cu60WHC_c1 <- "#1B9E77"

depression0Cu90WHC_c1 <- "#D95F02"

depression260Cu60WHC_c1 <- "#7570B3"

depression260Cu90WHC_c1 <-   "#2166AC"



##adding centriod
centroid_I_de <- cladeI_genus_nmds_points%>%
  subset(soil == "de")%>%
  group_by(treatment)%>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2) )




##NMDS plotting#####
####depression soil by treatment cladeI_genus_nmds_points with rel_abundance#####
nmds_I_de_reabun <- cladeI_genus_nmds_points%>%
  filter(soil == "de")%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment, fill = treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F) +
  geom_point() +
  coord_cartesian(xlim=c(-0.15, 0.05), ylim=c(0.05, 0.3)) +
  # geom_point(data = centroid_I_de, size = 5)+
  labs(title = "*nosZ* clade I",
    x="NMDS Axis 1",
    y="NMDS Axis 2",
    caption = "Depression soils, using relative abundance for metaMDS(k = 2)")+
  scale_color_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                              depression260Cu60WHC_c1,
                                              depression0Cu90WHC_c1,
                                              depression260Cu90WHC_c1),
                     labels=c("Depression soil Cu0 WHC60",
                              "Depression soil Cu260 WHC60",
                              "Depression soil Cu0 WHC90",
                              "Depression soil Cu260 WHC90"),guide = "none")+
  scale_fill_manual(name = NULL, values  = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                    guide = "none")+
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")



####depression soil by treatment cladeI_genus_nmds_points2 metaMDS(asv.bray1,  k = 2)#####
nmds_I_de_vegdist_reabun <- cladeI_genus_nmds_points2%>%
  filter(soil == "de")%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment, fill = treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F) +
  geom_point() +
  coord_cartesian(xlim=c(-0.06, 0), ylim=c(-0.02, 0.015)) +
  labs(title = "*nosZ* clade I",
    x="NMDS Axis 1",
    y="NMDS Axis 2",
    caption = "Depression soils, using output of vegdist (asv.bray1) for metaMDS(k = 2)") +
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
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.6, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")

nmds_I_de_reabun + nmds_I_de_vegdist_reabun

ggsave("~/Downloads/seq/MNDSplot/cladeI_de_vegdist_reabun_2.tiff",  width=8, height=4)

####cladeI_genus_nmds_points upland soil by treatment NMDS plot#####
cladeI_genus_nmds_points%>%
  filter(soil == "up")%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment, fill=treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.2, show.legend=F) +
  geom_point() +
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3.5, 2)) +
  labs(title = "*nosZ* clade I",
    x="NMDS Axis 1",
    y="NMDS Axis 2",
    caption = "Upland soils, using relative abundance for metaMDS(k = 2)")+
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
                                             upland260Cu90WHC_c1),
                    guide = "none")+
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")




####cladeI_genus_nmds_points, both soils by treatment#####
# Create a combined variable for treatment and soil
cladeI_genus_nmds_points <- cladeI_genus_nmds_points %>%
  mutate(treatment_soil = paste(treatment, soil, sep = "_"))%>%
  mutate(treatment_soil, factor(treatment_soil,
                                levels =c("upll_up",
                                          "uphl_up",
                                          "uplh_up",
                                          "uphh_up",
                                          "dell_de",
                                          "dehl_de",
                                          "delh_de",
                                          "dehh_de")))

# Define a custom color palette
color_palette <- c(
  "upll_up" = upland0Cu60WHC_c1,
  "uphl_up" = upland260Cu60WHC_c1,
  "uplh_up" = upland0Cu90WHC_c1,
  "uphh_up" = upland260Cu90WHC_c1,
  "dell_de" = depression0Cu60WHC_c1,
  "dehl_de" = depression260Cu60WHC_c1,
  "delh_de" = depression0Cu90WHC_c1,
  "dehh_de" = depression260Cu90WHC_c1)

# Define shape values
shape_values <- c(
  "upll_up" = 17,
  "uphl_up" = 17,
  "uplh_up" = 17,
  "uphh_up" = 17,
  "dell_de" = 19,
  "dehl_de" = 19,
  "delh_de" = 19,
  "dehh_de" = 19)

# Plot
ggplot(cladeI_genus_nmds_points,
       aes(x = NMDS1, y = NMDS2, color = treatment_soil, 
           fill = treatment_soil, shape = treatment_soil)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.75, alpha = 0.2, show.legend = F) +
  geom_point() +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3.6, 2)) +
  labs(title = "*nosZ* clade I",
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    caption = "All treatment, using relative abundance for metaMDS(k = 2)")+
  scale_color_manual(
    name = NULL,
    values = color_palette,
    labels = c(
      "upll_up" = "Upland Cu0 WHC60",
      "uphl_up" = "Upland Cu260 WHC60",
      "uplh_up" = "Upland Cu0 WHC90",
      "uphh_up" = "Upland Cu260 WHC90",
      "dell_de" = "Depression Cu0 WHC60",
      "dehl_de" = "Depression Cu260 WHC60",
      "delh_de" = "Depression Cu0 WHC90",
      "dehh_de" = "Depression Cu260 WHC90")) +
  scale_fill_manual(
    name = NULL,
    values = color_palette,
    labels = c(
      "upll_up" = "Upland Cu0 WHC60",
      "uphl_up" = "Upland Cu260 WHC60",
      "uplh_up" = "Upland Cu0 WHC90",
      "uphh_up" = "Upland Cu260 WHC90",
      "dell_de" = "Depression Cu0 WHC60",
      "dehl_de" = "Depression Cu260 WHC60",
      "delh_de" = "Depression Cu0 WHC90",
      "dehh_de" = "Depression Cu260 WHC90"),  guide = "none") +
  scale_shape_manual(
    name = NULL,
    values = shape_values,
    labels = c(
      "upll_up" = "Upland Cu0 WHC60",
      "uphl_up" = "Upland Cu260 WHC60",
      "uplh_up" = "Upland Cu0 WHC90",
      "uphh_up" = "Upland Cu260 WHC90",
      "dell_de" = "Depression Cu0 WHC60",
      "dehl_de" = "Depression Cu260 WHC60",
      "delh_de" = "Depression Cu0 WHC90",
      "dehh_de" = "Depression Cu260 WHC90")) +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.85, 0.2),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = -2, r = 3, b = 3, l = 3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b = 1, unit = "lines")),
    plot.margin = margin(l = 0.5, r = 0.5, t = 0.5, b = 0.5, unit = "lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")






##Clade II NMDS analysis#####
cII <-  read.table("~/Downloads/seq/total-spruce-taxa_clade_II.txt",
                   header=TRUE, sep="\t")%>%
  dplyr::select(taxa,clade,value,sample, soil, Cu, waterlevel)%>%
  mutate(treatment =  paste0(cII$soil,cII$Cu, cII$waterlevel))%>%
  filter(soil != "up") %>%
  group_by(sample)%>%
  mutate(rel_abund = (value / sum(value))*100)



meancII <-  cII %>%
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





# cladeII genus plot without NA new relative abundance

taxon_II_replicate  <- meancII %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 


taxon_poll_II_replicate <- taxon_II_replicate%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")



taxon_polled_II_replicate <- inner_join(taxon_II_replicate, taxon_poll_II_replicate, by = "taxon")%>%
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


tlist2 <- meancII_wide %>%
  select(sample, treatment, soil, Cu, waterlevel)%>%
  distinct(sample, treatment, soil, Cu, waterlevel)


ttaxon_polled_II_replicate <- taxon_polled_II_replicate%>%
  select(-mean)%>%
  pivot_wider(names_from = taxon, values_from = mean_rel_abund)%>%
  right_join( tlist2, by = "sample") 
#~asv.tax.t in the example
view(ttaxon_polled_II_replicate)


str(ttaxon_polled_II_replicate)
head(ttaxon_polled_II_replicate)



# Permanova test: are samples significantly different from each other 
# when grouped by site or by vegetation ?
library("vegan")
ttaxon_polled_II_replicate_imputed <- ttaxon_polled_II_replicate %>%
  replace_na(list(" Chitinophaga" = 0, " Cytophagaceae_bacterium"  = 0,
                  " Flavisolibacter" = 0, " Opitutaceae_bacterium" = 0,
                  " Pedobacter"   = 0,    " Pontibacter"  = 0,     " Prevotella" = 0, 
                  "Other"  = 0))





str(ttaxon_polled_II_replicate_imputed)

#colomn numebr
asv.bray1_II <- as.data.frame(as.matrix(vegdist(ttaxon_polled_II_replicate_imputed[,2:(ncol(ttaxon_polled_II_replicate_imputed)-4)], 
                                                diag = TRUE, upper = TRUE)))


adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$Cu, asv.bray1_II, permutations = 3000) # p = 0.2686
adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$waterlevel, asv.bray1_II, permutations = 3000) # p = 0.05265 .
adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$treatment, asv.bray1_II, permutations = 3000) # p =  0.0983 .

adonis2(asv.bray1_II ~ ttaxon_polled_II_replicate_imputed$Cu, asv.bray1_II, permutations = 3000, 
        strata = ttaxon_polled_II_replicate_imputed$waterlevel) # p = 0.2089

#####Run the NMDS####
# cladeII_genus_nmds <- ttaxon_polled_II_replicate_imputed%>%
#   column_to_rownames(var = "sample")%>%
#   metaMDS(ttaxon_polled_II_replicate_imputed[,2:(ncol(ttaxon_polled_II_replicate_imputed)-4)], k=2)
# 
# 
# cladeI_genus_nmds
# 
# summary(cladeI_genus_nmds)
# str(cladeI_genus_nmds)
# 
# 

# Create a new column for row names
ttaxon_polled_II_replicate_imputed$sample_for_row_names <- ttaxon_polled_II_replicate_imputed$sample

# Set the new column as row names
ttaxon_polled_II_replicate_imputed <- ttaxon_polled_II_replicate_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_II_replicate_imputed)

# Select only numeric columns for metaMDS
numeric_columns_II <- sapply(ttaxon_polled_II_replicate_imputed, is.numeric)
ttaxon_numeric_II <- ttaxon_polled_II_replicate_imputed[, numeric_columns_II]

# Run metaMDS on the numeric columns
cladeII_genus_nmds <- metaMDS(ttaxon_numeric_II, k = 2)

# View the results
print(cladeII_genus_nmds)


str(cladeII_genus_nmds)


# Stress plot
stressplot(cladeII_genus_nmds)

# svg("plots/asv_nmds_stressplot.svg", width = 4, height = 4)
stressplot(cladeII_genus_nmds)
dev.off()



###NMDS with asv.bray1_II###
cladeII_genus_nmds2 <- metaMDS(asv.bray1_II,  k = 2)
# Warning message:
#   In metaMDS(asv.bray1_II, k = 2) :
#   stress is (nearly) zero: you may have insufficient data

stressplot(cladeII_genus_nmds2)

str(cladeII_genus_nmds2)
print(cladeII_genus_nmds2)


#geting the points for plotting#
cladeII_genus_nmds_points <- as.data.frame(cladeII_genus_nmds$points)%>%
  mutate(sample = row.names(.))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist2, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))



cladeII_genus_nmds_points2 <- as.data.frame(cladeII_genus_nmds2$points)%>%
  mutate(sample = row.names(cladeII_genus_nmds$points))%>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2)%>%
  right_join(tlist2, by = "sample")%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))



str(cladeII_genus_nmds_points)


str(cladeII_genus_nmds_points2)


#####NMDS plotting#####
######depression soil by treatment cladeII_genus_nmds_points with rel_abundance#####
cladeII_genus_nmds_points%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment, fill = treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F) +
  geom_point() +
  coord_cartesian(xlim=c(-0.75, 1.05), ylim=c(-0.35, 0.35)) +
  labs(title = "*nosZ* clade II",
       x="NMDS Axis 1",
       y="NMDS Axis 2",
       caption = "Depression soils, using relative abundance for metaMDS(k = 2)")+
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
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot" )



######depression soil by treatment cladeII_genus_nmds_points2 metaMDS(asv.bray1,  k = 2)#####
cladeII_genus_nmds_points2%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment, fill = treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F) +
  geom_point() +
  # geom_richtext(data=my_legend,
  #               aes(x=x, y=y, label=label, color=color),
  #               hjust=0, lineheight=0.75,
  #               inherit.aes = FALSE, show.legend=FALSE,
  #               fill = NA, label.color = NA
  # ) +
  coord_cartesian(xlim=c(-0.6, 1.2), ylim=c(-0.0003, 0.0003)) +
  labs(title = "*nosZ* clade II",
       x="NMDS Axis 1",
       y="NMDS Axis 2",
       caption = "Depression soils, using output of vegdist (asv.bray1_II) for metaMDS(k = 2)") +
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
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.8, 0.15),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")



##Relative abundance genus stack bar plot####
#####cladeI ####
taxon_rel_abundI_genus2 <- meancI %>%
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

library(RColorBrewer)

colors_set3 <- brewer.pal(12, "Set3")
colors_paired <- brewer.pal(1, "Paired")  # Adding one color from another palette

# Combine the colors to have a total of 13
all_colors <- c(colors_set3, colors_paired)

inner_join(taxon_rel_abundI_genus2, taxon_pollIg2, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(treatment, taxon)%>%
  summarise(mean_rel_abund = sum(mean_rel_abund3),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1)) %>% 
  ggplot(aes(x=treatment, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL, values = all_colors) +
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
  scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt"))  



#######cladeII#####

meancII <-  cII %>%
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






taxon_pollIIg2 <- taxon_rel_abundII_genus2%>% 
  group_by(taxon)%>%
  summarise(pool = max(mean_rel_abund3) <3,
            mean = mean(mean_rel_abund3),
            .groups = "drop")

# Combine the colors to have a total of 13
library(RColorBrewer)
colors_set3 <- brewer.pal(12, "Set3")
colors_paired <- brewer.pal(1, "Paired")  # Adding one color from another palette
all_colors <- c(colors_set3, colors_paired)


#plotting
inner_join(taxon_rel_abundII_genus2, taxon_pollIIg2, by = "taxon")%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  group_by(treatment, taxon)%>%
  summarise(mean_rel_abund = sum(mean_rel_abund3),
            mean = min(mean),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon,mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1)) %>% 
  ggplot(aes(x=treatment, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL, values = all_colors) +
  scale_x_discrete(breaks=c("dehh",
                            "dehl",
                            "delh",
                            "dell"),
                   labels = c("depression<br>Cu260<br>water90",
                              "depression<br>Cu260<br>water60",
                              "depression<br>Cu0<br>water90",
                              "depression<br>Cu0<br>water60"
                   ))+
  scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt")) 



###cca####

tI <- ttaxon_polled_I_replicate_imputed[,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)]%>%
  as_tibble(rownames = "sample")

tI_matrix <- ttaxon_polled_I_replicate_imputed[,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)]

detI_matrix <- tI %>%
  inner_join(., tlist, by = "sample")%>%
  subset(soil == "de")%>%
  select(-sample,-treatment, -soil, -Cu, -waterlevel)%>%
  as.matrix()



###cca cladeI all ####

ord_cca_cI <- cca(tI_matrix~ soil + Cu, data = clade_env)
ord_cca_cIall <- cca(tI_matrix~ soil + Cu +waterlevel, data = clade_env)

plot(ord_cca_cI)
anova(ord_cca_cI)
# Model: cca(formula = tI_matrix ~ soil + Cu, data = clade_env)
# Df ChiSquare     F Pr(>F)    
# Model     2   0.26202 1.489  0.001 ***
#   Residual 26   2.28769              


plot(ord_cca_cIall)
anova(ord_cca_cIall)
# Model: cca(formula = tI_matrix ~ soil + Cu + waterlevel, data = clade_env)
# Df ChiSquare      F Pr(>F)   
# Model     3   0.35178 1.3337  0.004 **
#   Residual 25   2.19793       



#analyse terms separately
anova(ord_cca_cI, by="term", permutations=199)
anova(ord_cca_cIall, by="term", permutations=199)

#by marginal effects (“Type IIIeffects”):
anova(ord_cca_cI, by="mar", permutations=199)

#analyse significance of each axis
anova(ord_cca_cI, by="axis", permutations=199)

#envfit
envfit_I <- envfit(ord_cca_cIall, clade_env[,4:5], permu = 999, na.rm = TRUE)


##only depression soils cladeI
ord_cca_cIde <- cca(detI_matrix~ waterlevel + Cu, data = clade_envde)

plot(ord_cca_cIde)
anova(ord_cca_cIde)
# Model: cca(formula = detI ~ waterlevel + Cu, data = clade_envde)
# Df ChiSquare      F Pr(>F)  
# Model     2 0.0066872 2.0356  0.032 *
#   Residual 13 0.0213529 



#analyse terms separately
anova(ord_cca_cIde, by="term", permutations=199)


#by marginal effects (“Type IIIeffects”):
anova(ord_cca_cIde, by="mar", permutations=199)

#analyse significance of each axis
anova(ord_cca_cIde, by="axis", permutations=199)





###clade II cca####
tII_matrix <- ttaxon_polled_II_replicate_imputed[,2:(ncol(ttaxon_polled_II_replicate_imputed)-4)]

ord_cca_cII <- cca(tII_matrix~  waterlevel+Cu, data = clade_envde)

heatmap(dist_matrix_II)

plot(ord_cca_cII)
anova(ord_cca_cII)
# Model: cca(formula = tII_matrix ~ waterlevel + Cu, data = clade_envde)
# Df ChiSquare      F Pr(>F)  
# Model     2  0.097559 2.5956  0.034 *
#   Residual 13  0.244313                





#analyse terms separately
anova(ord_cca_cII, by="term", permutations=199)


#by marginal effects (“Type IIIeffects”):
anova(ord_cca_cII, by="mar", permutations=199)

#analyse significance of each axis
anova(ord_cca_cII, by="axis", permutations=199)

dist_matrix_I <- as.matrix(vegdist(ttaxon_polled_I_replicate_imputed[,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)], 
                                   diag = TRUE, upper = TRUE))


envfit_II <- envfit(ord_cca_cII2, clade_envde[,4:5], permu = 999, na.rm = TRUE)






##NMDS with only depression soil cladeI_genus_nmds_points_de with rel_abundance#####

# Create a new column for row names
ttaxon_polled_I_replicate_imputed$sample_for_row_names <- ttaxon_polled_I_replicate_imputed$sample

# Set the new column as row names
ttaxon_polled_I_replicate_imputed_de <- ttaxon_polled_I_replicate_imputed %>%
  filter(soil == "de")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed_de)

# Select only numeric columns for metaMDS
numeric_columns_de <- sapply(ttaxon_polled_I_replicate_imputed_de, is.numeric)
ttaxon_numeric_de <- ttaxon_polled_I_replicate_imputed_de[, numeric_columns_de]

# Run metaMDS on the numeric columns
cladeI_genus_nmds_de <- metaMDS(ttaxon_numeric_de, k = 2)


cladeI_genus_nmds_points_de <- as.data.frame(cladeI_genus_nmds_de$points)%>%
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





#NMD plotting
cladeI_genus_nmds_points_de%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2, color=treatment , fill = treatment)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F) +
  geom_point() +
  coord_cartesian(xlim=c(0.6, -0.8), ylim=c(0.6, -0.4)) +
  labs(title = "*nosZ* clade I",
       x="NMDS Axis 1",
       y="NMDS Axis 2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2)")+
  
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


####clade I depression envfit with cladeI_genus_nmds_de####
envfit_I_de_nmds <- envfit(cladeI_genus_nmds_de, clade_env1_de, permu = 999, 
                           na.rm = TRUE)


plot(cladeI_genus_nmds_de, display = "sites")
plot(envfit_I_de_nmds)



envfit_I_de_arrows_nmds <- as.data.frame(envfit_I_de_nmds$vectors$arrows) %>%
  mutate(varNames = rownames(.))
envfit_I_de_arrows_nmds <- envfit_I_de_arrows_nmds%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))



####significnat test nmds envfit cladeI de####
r_I_de_nmds_env <- as.matrix(envfit_I_de_nmds$vectors$r) 
p_I_de_nmds_env <- as.matrix(envfit_I_de_nmds$vectors$pvals)
env_p_I_de_nmds <- cbind(r_I_de_nmds_env,p_I_de_nmds_env)
colnames(env_p_I_de_nmds) <- c("r2","p_value")
K_I_de_nmds <- as.data.frame(env_p_I_de_nmds)%>%
  mutate(p_adjust = p.adjust(K_I_de_nmds$p_value, method = "BH"))
K_I_de_nmds

write_xlsx(K_I_de_nmds,"~/Downloads/seq/stat/K_I_de_nmds.xlsx")


####envfit cladeI_genus_nmds_de, and plotting######
#species
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





##NMDs species plotting + envfit 
cladeI_genus_nmds_points_de_sp%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2)) +
  geom_segment(data = cladeI_genus_nmds_points_de_sp1,
               aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black",
               alpha = 0.7) +
  geom_segment(data = envfit_I_de_arrows_nmds,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color= "#999999") +
  geom_label_repel(data = envfit_I_de_arrows_nmds,
                   mapping = aes(x = NMDS1, y = NMDS2, 
                                 label = varNames),
                   color="#0072b2",
                   min.segment.length = 0, seed = 42,
                   label.padding = unit(0, "lines"), 
                   label.size = NA,
                   parse = TRUE ) +
  geom_label_repel(data = cladeI_genus_nmds_points_de_sp1,
                   mapping = aes(x = NMDS1, y = NMDS2, label = sprintf("italic('%s')", taxon)),#label = paste0("italic(", taxon, ")" 
                   color="black",
                   alpha = 0.7,parse = TRUE
  ) +
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de, permu = 999, na.rm = TRUE)") +
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



##NMDs sites plotting +envfit 
cladeI_genus_nmds_points_de_si%>%
  ggplot(
    aes(x=NMDS1, y=NMDS2)) +
geom_segment(data = cladeI_genus_nmds_points_de_si1,
             aes(x = x_start, y = y_start, xend = NMDS1, yend = NMDS2),
             arrow = arrow(length = unit(0.2, "cm")), 
             color = "black",
             alpha = 0.7) +
  geom_segment(data = envfit_I_de_arrows_nmds,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color= "#999999") +
  geom_label_repel(data = envfit_I_de_arrows_nmds,
                   mapping = aes(x = NMDS1, y = NMDS2, label = varNames), 
                   color="#0072b2",
                   min.segment.length = 0, seed = 42, parse = TRUE
  ) +
  geom_label_repel(data = cladeI_genus_nmds_points_de_si1,
                   mapping = aes(x = NMDS1, y = NMDS2, label = taxon), 
                   color="black",
                   alpha = 0.7,parse = TRUE
  ) +
  labs(title = "*nosZ* clade I",
       x="NMDS1",
       y="NMDS2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de, permu = 999, na.rm = TRUE)") +
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


#######*NMDS ellipse + envfit clade I de####
cladeI_genus_nmds_points_de %>%
  ggplot(
    aes(x = NMDS1, y = NMDS2, color = treatment)) +
  stat_ellipse(geom = "polygon", aes(fill = treatment), 
               type = "norm", level = 0.75, alpha = 0.1, show.legend = FALSE) +
  geom_point() +
  coord_cartesian(xlim = c(-1,  1), ylim = c(-1.5,  0.5)) +
  labs(title = "*nosZ* clade I",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de, permu = 999, na.rm = TRUE)") +
  geom_segment(data = envfit_I_de_arrows_nmds,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color= "#999999") +
  geom_label_repel(data = envfit_I_de_arrows_nmds,
                   mapping = aes(x = NMDS1, y = NMDS2, label = varNames), 
                   color="#0072b2",
                   min.segment.length = 0, seed = 42, parse = TRUE
  ) +
  scale_color_manual(name = NULL, values = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                     labels = c("Depression soil Cu0 WHC60",
                                "Depression soil Cu260 WHC60",
                                "Depression soil Cu0 WHC90",
                                "Depression soil Cu260 WHC90")) +
  scale_fill_manual(name = NULL, values = c(depression0Cu60WHC_c1,
                                            depression260Cu60WHC_c1,
                                            depression0Cu90WHC_c1,
                                            depression260Cu90WHC_c1),
                    guide = "none") +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.75, 0.1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = -2, r = 3, b = 3, l = 3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b = 1, unit = "lines")),
    plot.margin = margin(l = 0.5, r = 0.5, t = 0.5, b = 0.5, unit = "lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")



###NMDs with only upland soils cladeI   #####

# Set the new column as row names
ttaxon_polled_I_replicate_imputed_up <- ttaxon_polled_I_replicate_imputed %>%
  filter(soil == "up")%>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_I_replicate_imputed_up)

# Select only numeric columns for metaMDS
numeric_columns_up <- sapply(ttaxon_polled_I_replicate_imputed_up, is.numeric)
ttaxon_numeric_up <- ttaxon_polled_I_replicate_imputed_up[, numeric_columns_up]



# Run metaMDS on the numeric columns
cladeI_genus_nmds_up <- metaMDS(ttaxon_numeric_up, k = 2)


####envfit cladeI_genus_nmds_up, and plotting######

envfit_I_up_nmds <- envfit(cladeI_genus_nmds_up, clade_env1_up, permu = 999, 
                           na.rm = TRUE)


plot(cladeI_genus_nmds_up, display = "sites")
plot(envfit_I_up_nmds)



####significnat test nmds envfit cladeI up####
r_I_up_nmds_env <- as.matrix(envfit_I_up_nmds$vectors$r) 
p_I_up_nmds_env <- as.matrix(envfit_I_up_nmds$vectors$pvals)
env_p_I_up_nmds <- cbind(r_I_up_nmds_env,p_I_up_nmds_env)
colnames(env_p_I_up_nmds) <- c("r2","p_value")
K_I_up_nmds <- as.data.frame(env_p_I_up_nmds)%>%
  mutate(p_adjust = p.adjust(K_I_up_nmds$p_value, method = "BH"))
K_I_up_nmds



##getting vectors for plotting
envfit_I_up_arrows_nmds <- as.data.frame(envfit_I_up_nmds$vectors$arrows) %>%
  mutate(varNames = rownames(.))

envfit_I_up_arrows_nmds <- envfit_I_up_arrows_nmds%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      TRUE ~ varNames
    ))




#species
cladeI_genus_nmds_points_up_sp <- as.data.frame(scores(cladeI_genus_nmds_up)$species)%>%
  mutate(sample = row.names(.))



cladeI_genus_nmds_points_up_sp1 <- cladeI_genus_nmds_points_up_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "(.*)_group", "*\\1* group"))
  


#site
cladeI_genus_nmds_points_up_si <- as.data.frame(scores(cladeI_genus_nmds_up)$site)%>%
  mutate(sample = row.names(.))


cladeI_genus_nmds_points_up_si1 <- cladeI_genus_nmds_points_up_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")



#######*NMDS ellipse + envfit clade I up####
cladeI_genus_nmds_points_up <- as.data.frame(cladeI_genus_nmds_up$points)%>%
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



cladeI_genus_nmds_points_up %>%
  ggplot(
    aes(x = NMDS1, y = NMDS2, color = treatment)) +
  stat_ellipse(geom = "polygon", aes(fill = treatment), type = "norm", level = 0.75, alpha = 0.1, show.legend = FALSE) +
  geom_point() +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 1.5)) +
  labs(title = "*nosZ* clade I",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2",
       caption = "Upland soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_up, clade_env1_up, permu = 999, na.rm = TRUE)") +
  geom_segment(data = envfit_I_up_arrows_nmds,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color= "#999999") +
  geom_label_repel(data = envfit_I_up_arrows_nmds,
                   mapping = aes(x = NMDS1, y = NMDS2, label = varNames), 
                   color="#0072b2",
                   min.segment.length = 0, seed = 42, parse = TRUE
  ) +
  scale_color_manual(name = NULL, values = c(upland0Cu60WHC_c1,
                                             upland260Cu60WHC_c1,
                                             upland0Cu90WHC_c1,
                                             upland260Cu90WHC_c1),
                     labels = c("upland soil Cu0 WHC60",
                                "upland soil Cu260 WHC60",
                                "upland soil Cu0 WHC90",
                                "upland soil Cu260 WHC90")) +
  scale_fill_manual(name = NULL, values = c(upland0Cu60WHC_c1,
                                            upland260Cu60WHC_c1,
                                            upland0Cu90WHC_c1,
                                            upland260Cu90WHC_c1),
                    guide = "none") +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.75, 0.15),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = -2, r = 3, b = 3, l = 3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b = 1, unit = "lines")),
    plot.margin = margin(l = 0.5, r = 0.5, t = 0.5, b = 0.5, unit = "lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")


##clade II ####
ttaxon_polled_II_replicate_imputed$sample_for_row_names <- ttaxon_polled_II_replicate_imputed$sample

# Set the new column as row names
ttaxon_polled_II_replicate_imputed <- ttaxon_polled_II_replicate_imputed %>%
  column_to_rownames(var = "sample_for_row_names")

# Check the structure to identify non-numeric columns
str(ttaxon_polled_II_replicate_imputed)

# Select only numeric columns for metaMDS
numeric_columns_II <- sapply(ttaxon_polled_II_replicate_imputed, is.numeric)
ttaxon_numeric_II <- ttaxon_polled_II_replicate_imputed[, numeric_columns_II]

# Run metaMDS on the numeric columns
cladeII_genus_nmds <- metaMDS(ttaxon_numeric_II, k = 2)


envfit_II_de_nmds <- envfit(cladeII_genus_nmds, clade_env1_de, permu = 999, 
                            na.rm = TRUE)


plot(cladeII_genus_nmds, display = "sites")
plot(envfit_II_de_nmds)


##significnat test nmds envfit cladeI de####
r_II_de_nmds_env <- as.matrix(envfit_II_de_nmds$vectors$r) 
p_II_de_nmds_env <- as.matrix(envfit_II_de_nmds$vectors$pvals)
env_p_II_de_nmds <- cbind(r_II_de_nmds_env,p_II_de_nmds_env)
colnames(env_p_II_de_nmds) <- c("r2","p_value")
K_II_de_nmds <- as.data.frame(env_p_II_de_nmds)%>%
  mutate(p_adjust = p.adjust(K_II_de_nmds$p_value, method = "BH"))
K_II_de_nmds




envfit_II_de_arrows_nmds <- as.data.frame(envfit_II_de_nmds$vectors$arrows) %>%
  mutate(varNames = rownames(.))
envfit_II_de_arrows_nmds <- envfit_II_de_arrows_nmds%>%
  mutate(
    varNames = case_when(
      varNames == "NO3" ~ as.character(parse(text = "NO[3]^\"-\"")),
      varNames == "NH4" ~ as.character(parse(text = "NH[4]^\"+\"")),
      varNames == "Cumass" ~ as.character(c("Cu")),
      
      TRUE ~ varNames))




##envfit cladeII_genus_nmds_de, and plotting######
#sp
cladeII_genus_nmds_points_de_sp <- as.data.frame(scores(cladeII_genus_nmds)$species)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_sp1 <- cladeII_genus_nmds_points_de_sp %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")%>%
  mutate(taxon = str_replace(taxon,
                             "unclassified_(.*)", "Unclassified \\1"),
         taxon = str_replace(taxon,
                             "(.*)_group", "\\1 ")) 



#site
cladeII_genus_nmds_points_de_si <- as.data.frame(scores(cladeII_genus_nmds)$site)%>%
  mutate(sample = row.names(.))


cladeII_genus_nmds_points_de_si1 <- cladeII_genus_nmds_points_de_si %>%
  mutate(x_start = 0, y_start = 0)%>%
  as_tibble(rownames="taxon")  




cladeII_genus_nmds_points %>%
  ggplot(
    aes(x = NMDS1, y = NMDS2, color = treatment)) +
  stat_ellipse(geom = "polygon", aes(fill = treatment), 
               type = "norm", level = 0.75, alpha = 0.1, show.legend = FALSE) +
  geom_point() +
  coord_cartesian(xlim = c(-0.8,  1), ylim = c(-1.5,  1)) +
  labs(title = "*nosZ* clade II",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2",
       caption = "Depression soils, using depression soil data relative abundance for metaMDS(k = 2),
       envfit(cladeI_genus_nmds_de, clade_env1_de, permu = 999, na.rm = TRUE)") +
  geom_segment(data = envfit_II_de_arrows_nmds,
               mapping = aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color= "#999999") +
  geom_label_repel(data = envfit_II_de_arrows_nmds,
                   mapping = aes(x = NMDS1, y = NMDS2, label = varNames), 
                   color="#0072b2",
                   min.segment.length = 0, seed = 42, parse = TRUE
  ) +
  scale_color_manual(name = NULL, values = c(depression0Cu60WHC_c1,
                                             depression260Cu60WHC_c1,
                                             depression0Cu90WHC_c1,
                                             depression260Cu90WHC_c1),
                     labels = c("Depression soil Cu0 WHC60",
                                "Depression soil Cu260 WHC60",
                                "Depression soil Cu0 WHC90",
                                "Depression soil Cu260 WHC90")) +
  scale_fill_manual(name = NULL, values = c(depression0Cu60WHC_c1,
                                            depression260Cu60WHC_c1,
                                            depression0Cu90WHC_c1,
                                            depression260Cu90WHC_c1),
                    guide = "none") +
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(0.75, 0.1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = -2, r = 3, b = 3, l = 3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b = 1, unit = "lines")),
    plot.margin = margin(l = 0.5, r = 0.5, t = 0.5, b = 0.5, unit = "lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot")
