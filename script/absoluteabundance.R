###clade I#####


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
cIcp <- cI%>%
  left_join(cpclade, by = "sample")

# absolute abundance of nosZ community members, aab
meancIcp <-  cIcp %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment, cladeI)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100,
         N = N)%>%
  mutate(aab = rel_abund *cladeI/100)%>%
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund,N, treatment, cladeI, aab)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
               names_to="level",
               values_to="taxon")


taxon_I_replicate_cp  <- meancIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 


hist(taxon_I_replicate_cp$aab)

taxon_poll_I_replicate_cp <- taxon_I_replicate_cp%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")




taxon_aabI_genus0 <- meancIcp %>%
  # mutate(aab = rel_abund *cladeI)%>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(aab_sum = sum(aab),
            mean_aab = mean(aab_sum), .groups="drop")%>%
  inner_join( taxon_poll_I_replicate_cp, by = "taxon")


taxon_aabI_genus <- taxon_aabI_genus0%>%
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
  mutate(mean_aab2 = mean_aab/replicates)%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))


# taxon_poll_I_replicate_cp


# inner_join(taxon_rel_abundI_family, taxon_poll_I_replicate_cp, by = "taxon")%>%

#Neisseria should also be pooled 
taxon_aabI_genus <- taxon_aabI_genus%>%
  mutate(pool = if_else(taxon == "*Neisseria*", TRUE, pool),
         pool = if_else(taxon == "*Polaromonas*", TRUE, pool))



taxon_aabI_genus_plot <-   taxon_aabI_genus%>%
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
  summarise(mean_rel_aab = sum(mean_aab2),
            mean = min(mean_aab2),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE), 
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1),
         soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"),
         soil = factor(soil, levels = c("u", "d")),
         treatment1 = str_replace(treatment,"^\\w*(\\w)(\\w)", "\\1\\2"),
         treatment1 = factor(treatment1,
                             levels = 
                               c("ll",
                                 "hl",
                                 "lh",
                                 "hh"
                               )))


#####plotting #####
all_colors2 <- c( "#A6CEE3", "#FFFFB3",   "#FB8072", "#7570B3",
                  "#FDB462" ,"#B3DE69" ,"#FCCDE5", "#BC80BD" , "#FFED6F", "#E5C494","gray","#8DD3C7", "#1F78B4","#FED9A6")

all_colors2cp <- c( "#A6CEE3", "#FFFFB3" ,"#B3DE69","#FCCDE5","#7570B3",
                    "#FDB462", "#FFED6F","#FB8072",  "#BC80BD" , "#E5C494" ,"gray", "#1F78B4", "#8DD3C7")


all_colors3cp <- c( "#A6CEE3", "#FFFFB3" ,"#B3DE69","#FCCDE5","#7570B3",
                    "#FDB462", "#FFED6F","#FB8072",  "#BC80BD" , "#E5C494","#FED9A6" ,"gray", "#1F78B4", "#8DD3C7")


all_colors3 <- c( "#A6CEE3", "#FFFFB3",   "#FB8072", "#7570B3",
                  "#FDB462" ,"#B3DE69" ,"#FCCDE5", "#BC80BD" , "#FFED6F", "#E5C494","gray","#8DD3C7","#BEBADA", "#1F78B4","#FED9A6")


taxon_aabI_genus_plot %>%
  filter(soil  == "u")%>%
  ggplot(aes(x=treatment1, y=mean_rel_aab, fill=taxon)) +
  # facet_grid(.~factor(soil, labels = c("u" = "upland soil", 
  #                                      "d" = "depression soil")))+
  geom_col() +
  scale_fill_manual(name=NULL,
                    # breaks=c("*Afipia*", "*Azospirillum*","*Bradyrhizobium*","*Hyphomicrobium*",
                    #          "*Mesorhizobium*", "*Microvirga*", "*Noviherbaspirillum*",
                    #          "*Pseudomonas*","*Ramlibacter*","*Rhizobium/Agrobacterium* group",
                    #          "*Skermanella*", "Unclassified *Acetobacteraceae*",
                    #          "Other"),
                    values = all_colors3cp) +
  #the last color "gray" --> "Other"
  # scale_fill_manual(name = NULL, 
  #                   
  #                   values = all_colors) +
  scale_y_continuous(labels = scales::scientific_format())+
  scale_x_discrete(breaks=c(
    "ll",
    "hl",
    "lh",
    "hh"),
    labels = c("60% WHC<br>Cu 0 mM",
               "60% WHC<br>Cu 260 mM",
               "90% WHC<br>Cu 0 mM",
               "90% WHC<br>Cu 260 mM"))+
  # labels=c("Healthy",
  #          "Diarrhea,<br>*C. difficile*<br>negative",
  #          "Diarrhea,<br>*C. difficile*<br>positive")) +
  # scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Absolute abundance of *nosZ* cladeI community") +
  theme_classic() +
  theme(strip.text = element_markdown(size = 18),
        strip.background =element_rect(color = NA),
        axis.text.x = element_markdown(size = 12),
        axis.title.y  = element_markdown(size = 16),
        axis.text.y = element_markdown(size = 14),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(28, "pt"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

ggsave("~/Downloads/seq/plot/aab/cladeI_genusXNA_up.tiff",  width=10, height=8)
ggsave("~/Downloads/seq/plot/aab/cladeI_genusXNA_up.jpeg",  width=10, height=8)
ggsave("~/Downloads/seq/aab/plot/cladeI_genusXNA_up.jpeg",  width=10, height=8)


taxon_aabI_genus_plot %>%
  filter(soil  == "d")%>%
  ggplot(aes(x=treatment1, y=mean_rel_aab, fill=taxon)) +
  # facet_grid(.~factor(soil, labels = c("u" = "upland soil", 
  #                                      "d" = "depression soil")))+
  geom_col() +
  scale_fill_manual(name=NULL,
                    # breaks=c("*Afipia*", "*Azospirillum*","*Bradyrhizobium*","*Hyphomicrobium*",
                    #          "*Mesorhizobium*", "*Microvirga*", "*Noviherbaspirillum*",
                    #          "*Pseudomonas*","*Ramlibacter*","*Rhizobium/Agrobacterium* group",
                    #          "*Skermanella*", "Unclassified *Acetobacteraceae*",
                    #          "Other"),
                    values = all_colors2cp) +
  #the last color "gray" --> "Other"
  # scale_fill_manual(name = NULL, 
  #                   
  #                   values = all_colors) +
  scale_y_continuous(labels = scales::scientific_format())+
  scale_x_discrete(breaks=c(
    "ll",
    "hl",
    "lh",
    "hh"),
    labels = c("60% WHC<br>Cu 0 mM",
               "60% WHC<br>Cu 260 mM",
               "90% WHC<br>Cu 0 mM",
               "90% WHC<br>Cu 260 mM"))+
  # labels=c("Healthy",
  #          "Diarrhea,<br>*C. difficile*<br>negative",
  #          "Diarrhea,<br>*C. difficile*<br>positive")) +
  # scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Absolute abundance of *nosZ* cladeI community") +
  theme_classic() +
  theme(strip.text = element_markdown(size = 18),
        strip.background =element_rect(color = NA),
        axis.text.x = element_markdown(size = 12),
        axis.title.y  = element_markdown(size = 16),
        axis.text.y = element_markdown(size = 14),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(28, "pt"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

ggsave("~/Downloads/seq/plot/aab/cladeI_genusXNA_de.tiff",  width=10, height=8)
ggsave("~/Downloads/seq/plot/aab/cladeI_genusXNA_de.jpeg",  width=10, height=8)
ggsave("~/Downloads/seq/aab/plot/cladeI_genusXNA_de.jpeg",  width=10, height=8)


###clade II #####


cII <-  read.table("~/Downloads/seq/total-spruce-taxa_clade_II.txt",
                   header=TRUE, sep="\t")%>%
  dplyr::select(taxa,clade,value,sample, soil, Cu, waterlevel)

cII <- cII%>%
  mutate(treatment =  paste0(cII$soil,cII$Cu, cII$waterlevel))%>%
  filter(soil != "up") %>%
  group_by(sample)%>%
  mutate(rel_abund = (value / sum(value))*100)



cpcladeII <- Cu_gene%>%
  filter(soil == "depression")%>%
  mutate(sample = as.character(ID),
        cladeII = nosz_cladeII )%>%
  mutate(sample = str_c("z", sample))%>%
  dplyr::select(sample , cladeII)%>%
  mutate(cladeII = as.numeric(cladeII))




cIIcp <- cII%>%
  filter(sample != "z50")%>%
  filter(sample != "z72")%>%
  left_join(cpcladeII, by = "sample")


str(cIIcp)



meancIIcp <-  cIIcp %>%
  group_by(taxa, soil, Cu, waterlevel, sample, treatment,cladeII)%>%
  summarize(N=n(),
            .groups = "drop")%>%
  separate(taxa,
           into=c("root","kingdom", "phylum", "class", "order", "family", "genus","species"),
           sep=";")%>%
  group_by(sample)%>%
  mutate(rel_abund = (N / sum(N))*100,
         N = N)%>%
  mutate(aab = rel_abund *cladeII/100)%>%
  dplyr::select(kingdom, phylum, class, order, 
                family, genus, species, sample,soil, Cu, waterlevel, rel_abund,N, treatment, cladeII, aab)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
               names_to="level",
               values_to="taxon")  




taxon_II_replicate_cp  <- meancIIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 


hist(taxon_II_replicate_cp$aab)



taxon_poll_II_replicate_cp <- taxon_II_replicate_cp%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")



taxon_aabII_genus0 <- meancIIcp %>%
  filter(soil == "de")%>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  group_by(treatment, taxon) %>%
  summarize(aab_sum = sum(aab),
            mean_aab = mean(aab_sum), .groups="drop")%>%
  inner_join( taxon_poll_II_replicate_cp, by = "taxon")




taxon_aabII_genus <- taxon_aabII_genus0%>%
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
    treatment == "dehh" ~ 3,
    treatment == "delh" ~ 4,
    treatment == "dell" ~ 3,
    treatment == "dehl" ~ 4))%>%
  mutate(mean_aab2 = mean_aab/replicates)%>%
  mutate(soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"))


taxon_aabII_genus <- taxon_aabII_genus%>%
  mutate(pool = if_else(taxon == "*Prevotella*", TRUE, pool))


taxon_aabII_genus_plot <-   taxon_aabII_genus%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))%>%
  group_by(treatment, taxon)%>%
  summarise(mean_rel_aab = sum(mean_aab2),
            mean = min(mean_aab2),
            .groups = "drop")%>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE), 
         taxon = fct_shift(taxon, n = 1)) %>%
  #if turn off desc, Firmicutes will be on the bottom
  mutate(taxon = factor(taxon),
         taxon = fct_shift(taxon, n = 1),
         soil  = str_replace(treatment,"^(\\w)\\w*", "\\1"),
         soil = factor(soil, levels = c("u", "d")),
         treatment1 = str_replace(treatment,"^\\w*(\\w)(\\w)", "\\1\\2"),
         treatment1 = factor(treatment1,
                             levels = 
                               c("ll",
                                 "hl",
                                 "lh",
                                 "hh"
                               )))



all_colors_II2 <- c(  "#1F78B4", "#FDB462","#FCCDE5",
                      "#B3DE69","gray" , "#A6CEE3","#E5C494")


all_colors_II3 <- c(  "#1F78B4", "#FDB462","#FCCDE5",
                      "#B3DE69", "#BC80BD","gray" , "#A6CEE3","#E5C494")



taxon_aabII_genus_plot %>% 
  ggplot(aes(x=treatment, y=mean_rel_aab, fill=taxon)) +
  geom_col() +
  # scale_fill_manual(name=NULL,
  #                   breaks = c("*Alphaproteobacteria*", "*Betaproteobacteria*",
  #                              "*Gammaproteobacteria*",
  #                              "*Unclassified*Bacteria**"),
  #                   values  = c(brewer.pal(4, "Dark2"))) +
  scale_fill_manual(name = NULL, 
                    # breaks = c("*Chitinophaga*", "*Cytophagaceae bacterium*",
                    #            "*Flavisolibacter*", "*Opitutaceae bacterium*",
                    #            "*Pedobacter*", "*Pontibacter*", "Other"),
                    
                    values = all_colors_II2) +
  scale_x_discrete(breaks=c("dell",
                            "dehl",
                            "delh",
                            "dehh"),
                   labels = c(
                     "60% WHC<br>Cu 0 mM",
                     "60% WHC<br>Cu 260 mM",
                     "90% WHC<br>Cu 0 mM",
                     "90% WHC<br>Cu 260 mM"))+
  scale_y_continuous(labels = scales::scientific_format())+
  # scale_y_continuous(expand = c(0 , 0))+
  labs(x=NULL,
       y="Absolute abundance of *nosZ* cladeII community") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size = 12),
        axis.text.y = element_markdown(size = 12),
        legend.text = element_markdown(size = 12),
        legend.key.size = unit(28, "pt"),
        axis.title.y  = element_markdown(size = 14),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


ggsave("~/Downloads/seq/plot/aab/cladeII_genusXNA_de.tiff",  width=10, height=8)

ggsave("~/Downloads/seq/plot/aab/cladeII_genusXNA_de.jpeg",  width=10, height=8)
ggsave("~/Downloads/seq/aab/plot/cladeII_genusXNA_de.jpeg",  width=10, height=8)




###correlation with N2O#####
###clade I#####
####pivote data#####

taxon_I_replicate_cp_cr  <- meancIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 

taxon_poll_I_replicate_cp_cr <- taxon_I_replicate_cp%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")%>%
  mutate(pool = if_else(taxon == " Neisseria", TRUE, pool),
         pool = if_else(taxon == " Polaromonas", TRUE, pool))




taxon_aabI_genus0_cr <- meancIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  inner_join(taxon_poll_I_replicate_cp_cr, by = "taxon")




taxon_aabI_genus_cr0 <- taxon_aabI_genus0_cr%>%
  mutate(taxon = str_trim(taxon),
         taxon = str_replace_all(taxon, "\\s", ""))  # Remove any whitespace

Cu_N2O_cumuf_innerjoin <- Cu_N2O_cumuf %>%
  mutate(sample = str_c("n", ID))%>%
  select(flux, sample)

taxon_aabI_genus_cr <- taxon_aabI_genus_cr0%>%
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
  inner_join(Cu_N2O_cumuf_innerjoin, by = "sample")%>%
  group_by(sample, taxon, soil,  Cu,    waterlevel)%>%
  summarise(rel_abundsum = sum(rel_abund2), #rel_abund --> with NAs, rel_abund2 --> without NAs
            N = sum(N),
            aab = sum(aab))





taxon_aabI_genus_cr_up <- taxon_aabI_genus_cr%>%
  filter(soil =="up")

taxon_aabI_genus_cr_de <- taxon_aabI_genus_cr%>%
  filter(soil =="de")

str(taxon_aabI_genus_cr_up)

t <- taxon_aabI_genus_cr_up %>%
  filter(taxon == "Other")%>%
  pivot_wider(names_from = taxon, values_from = aab)



library(purrr)
# Create a named list of pivoted data for each unique taxon, into a list
pivoted_data_list <- unique(taxon_aabI_genus_cr_up$taxon) %>%
  set_names() %>%
  map(~ taxon_aabI_genus_cr_up %>%
        filter(taxon == .x) %>%
        pivot_wider(names_from = taxon, values_from = aab))



pivoted_data <- NULL
# Loop through eNULL# Loop through each unique value in taxon
for (i in unique(taxon_aabI_genus_cr_up$taxon)) {
  # Filter and pivot wider
  pivoted_data <- taxon_aabI_genus_cr_up %>%
    filter(taxon == i) %>%
    pivot_wider(names_from = taxon, values_from = aab)%>%
    inner_join(Cu_N2O_cumuf_innerjoin, by = "sample")
  
  # Dynamically assign the dataset to a variable
  assign(paste0("pivoted_I_up_", i), pivoted_data)
}




pivoted_data <- NULL
# Loop through eNULL# Loop through each unique value in taxon
for (i in unique(taxon_aabI_genus_cr_de$taxon)) {
  # Filter and pivot wider
  pivoted_data <- taxon_aabI_genus_cr_de %>%
    filter(taxon == i) %>%
    pivot_wider(names_from = taxon, values_from = aab)%>%
    inner_join(Cu_N2O_cumuf_innerjoin, by = "sample")
  
  # Dynamically assign the dataset to a variable
  assign(paste0("pivoted_I_de_", i), pivoted_data)
}





view(taxon_aabI_genus_cr_up)




####select correlation up data #####
datasets <- ls(pattern = "^pivoted_I_up_")

# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)%>%
    ungroup() %>%        # Remove grouping structure
    select(-sample, -Cu, -soil)
  
  # Explicitly ungroup the data to avoid adding missing grouping variables
  if ("group_vars" %in% names(attributes(data))) {
    data <- data %>% ungroup()
  }
  
  # Select the required columns (i and flux)
  selected_data <- data %>% 
    
    select(matches("i"), flux)  # Adjust column name selection logic as needed
  
  # Save the new dataset with "_cor" suffix
  assign(paste0(dataset, "_cor"), selected_data)
}


pivoted_I_up_Other_cor <- pivoted_I_up_Other%>%
  ungroup() %>%        # Remove grouping structure
  select(Other, flux)

pivoted_I_up_Pseudomonas_cor <- pivoted_I_up_Pseudomonas%>%
  ungroup() %>%        # Remove grouping structure
  select(Pseudomonas, flux)

pivoted_I_up_Skermanella_cor <- pivoted_I_up_Skermanella%>%
  ungroup() %>%        # Remove grouping structure
  select(Skermanella, flux)



####select correlation de data #####
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_I_de_")

# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)%>%
    ungroup() %>%        # Remove grouping structure
  select(-sample, -Cu, -soil)
  
  # Explicitly ungroup the data to avoid adding missing grouping variables
  if ("group_vars" %in% names(attributes(data))) {
    data <- data %>% ungroup()
  }
  
  # Select the required columns (i and flux)
  selected_data <- data %>% 
    
    select(matches("i"), flux)  # Adjust column name selection logic as needed
  
  # Save the new dataset with "_cor" suffix
  assign(paste0(dataset, "_cor"), selected_data)
}


pivoted_I_de_Other_cor <- pivoted_I_de_Other%>%
  ungroup() %>%        # Remove grouping structure
  select(Other, flux)

pivoted_I_de_Pseudomonas_cor <- pivoted_I_de_Pseudomonas%>%
  ungroup() %>%        # Remove grouping structure
  select(Pseudomonas, flux)

pivoted_I_de_Skermanella_cor <- pivoted_I_de_Skermanella%>%
  ungroup() %>%        # Remove grouping structure
  select(Skermanella, flux)




####correlation up cladeI####
##spearman
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_I_up_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_I_up <- NULL
t_p_I_up <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "spearman")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)
  
  t_r_I_up <- rbind(t_r_I_up, r_data)       
  t_p_I_up <- rbind(t_p_I_up, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_I_up_Other_cor"]]


t_r_I_up_df <- as.data.frame(t_r_I_up)%>%
  rownames_to_column(var = "row")
write_xlsx(t_r_I_up_df,"~/Downloads/seq/aab/cor/aab_I_up_r.xlsx")  


t_p_I_up_df <- as.data.frame(t_p_I_up)%>%
  rownames_to_column(var = "row")
write_xlsx(t_p_I_up_df,"~/Downloads/seq/aab/cor/aab_I_up_p.xlsx")  



##pearson
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_I_up_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_I_up_pear <- NULL
t_p_I_up_pear <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "pearson")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)
  
  t_r_I_up_pear <- rbind(t_r_I_up_pear, r_data)       
  t_p_I_up_pear <- rbind(t_p_I_up_pear, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_I_up_Other_cor"]]


t_r_I_up_df_pear <- as.data.frame(t_r_I_up_pear)%>%
  rownames_to_column(var = "row")

write_xlsx(t_r_I_up_df_pear,"~/Downloads/seq/aab/cor/aab_I_up_r_pear.xlsx")  


t_p_I_up_df_pear <- as.data.frame(t_p_I_up_pear)%>%
  rownames_to_column(var = "row")
write_xlsx(t_p_I_up_df_pear,"~/Downloads/seq/aab/cor/aab_I_up_p_pear.xlsx")  


pivoted_Afipia%>%
  subset(soil == "up")%>%
  ggplot(aes(x = Afipia, y = flux))+
  geom_point()+
  geom_smooth(method = "lm") 


pivoted_I_de_Afipia%>%
  ggplot(aes(x = Afipia, y = flux))+
  geom_point()+
  geom_smooth(method = "lm")  # Add a linear regression line

pivoted_I_de_Azospirillum%>%
  ggplot(aes(x = Azospirillum, y = flux))+
  geom_point()+
  geom_smooth(method = "lm")


pivoted_I_up_Azospirillum%>%
  ggplot(aes(x = Azospirillum, y = flux))+
  geom_point()+
  geom_smooth(method = "lm")


####correlation de clade I####
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_I_de_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_I_de <- NULL
t_p_I_de <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "spearman")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)

  t_r_I_de <- rbind(t_r_I_de, r_data)       
  t_p_I_de <- rbind(t_p_I_de, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_I_de_Other_cor"]]


t_r_I_de_df <- as.data.frame(t_r_I_de)%>%
  rownames_to_column(var = "row")
write_xlsx(t_r_I_de_df,"~/Downloads/seq/aab/cor/aab_I_de_r.xlsx")  


t_p_I_de_df <- as.data.frame(t_p_I_de)%>%
  rownames_to_column(var = "row")
write_xlsx(t_p_I_de_df,"~/Downloads/seq/aab/cor/aab_I_de_p.xlsx")  


###pearson
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_I_de_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_I_de_pear <- NULL
t_p_I_de_pear <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "pearson")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)
  
  t_r_I_de_pear <- rbind(t_r_I_de_pear, r_data)       
  t_p_I_de_pear <- rbind(t_p_I_de_pear, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_I_de_Other_cor"]]


t_r_I_de_df_pear <- as.data.frame(t_r_I_de_pear)%>%
  rownames_to_column(var = "row")
write_xlsx(t_r_I_de_df_pear,"~/Downloads/seq/aab/cor/aab_I_de_r_pear.xlsx")  


t_p_I_de_df_pear <- as.data.frame(t_p_I_de_pear)%>%
  rownames_to_column(var = "row")
write_xlsx(t_p_I_de_df_pear,"~/Downloads/seq/aab/cor/aab_I_de_p_pear.xlsx")  







###clade II####
####pivote data#####
taxon_II_replicate_cp_cr  <- meancIIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop") 


taxon_poll_II_replicate_cp_cr <- taxon_II_replicate_cp%>%
  group_by(taxon)%>%
  summarise(pool = max(rel_abund2) <3,
            mean = mean(rel_abund2),
            .groups = "drop")%>%
  mutate(pool = if_else(taxon == " Prevotella", TRUE, pool))




taxon_aabII_genus0_cr <- meancIIcp %>%
  filter(level=="genus") %>%
  filter(!is.na(taxon)) %>%
  group_by(sample) %>%
  mutate(rel_abund2 = (N / sum(N))*100, .groups ="drop")%>%
  inner_join(taxon_poll_II_replicate_cp_cr, by = "taxon")


taxon_aabII_genus_cr0 <- taxon_aabII_genus0_cr%>%
  mutate(taxon = str_trim(taxon),
         taxon = str_replace_all(taxon, "\\s", ""))  # Remove any whitespace

Cu_N2O_cumuf_innerjoin_II <- Cu_N2O_cumuf %>%
  mutate(sample = str_c("z", ID))%>%
  select(flux, sample)


taxon_aabII_genus_cr <- taxon_aabII_genus_cr0%>%
  mutate(taxon = if_else(pool, "Other", taxon))%>%
  mutate(treatment = factor(treatment,
                            levels = 
                              c(
                                "dell",
                                "dehl",
                                "delh",
                                "dehh")))%>%
  inner_join(Cu_N2O_cumuf_innerjoin_II, by = "sample")%>%
  group_by(sample, taxon, soil,  Cu,    waterlevel)%>%
  summarise(rel_abundsum = sum(rel_abund2), #rel_abund --> with NAs, rel_abund2 --> without NAs
            N = sum(N),
            aab = sum(aab))




pivoted_data <- NULL
# Loop through eNULL# Loop through each unique value in taxon
for (i in unique(taxon_aabII_genus_cr$taxon)) {
  # Filter and pivot wider
  pivoted_data <- taxon_aabII_genus_cr %>%
    filter(taxon == i) %>%
    pivot_wider(names_from = taxon, values_from = aab)%>%
    inner_join(Cu_N2O_cumuf_innerjoin_II, by = "sample")
  
  # Dynamically assign the dataset to a variable
  assign(paste0("pivoted_II_de_", i), pivoted_data)
}




####select correlation clade II de data #####
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_II_de_")

# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)%>%
    ungroup() %>%        # Remove grouping structure
    select(-sample, -Cu, -soil)
  
  # Explicitly ungroup the data to avoid adding missing grouping variables
  if ("group_vars" %in% names(attributes(data))) {
    data <- data %>% ungroup()
  }
  
  # Select the required columns (i and flux)
  selected_data <- data %>% 
    
    select(matches("i"), flux)  # Adjust column name selection logic as needed
  
  # Save the new dataset with "_cor" suffix
  assign(paste0(dataset, "_cor"), selected_data)
}


pivoted_II_de_Other_cor <- pivoted_II_de_Other%>%
  ungroup() %>%        # Remove grouping structure
  select(Other, flux)

pivoted_II_de_Pedobacter_cor <- pivoted_II_de_Pedobacter%>%
  ungroup() %>%        # Remove grouping structure
  select(Pedobacter, flux)



####correlation de clade II####
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_II_de_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_II_de <- NULL
t_p_II_de <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "spearman")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)
  
  t_r_II_de <- rbind(t_r_II_de, r_data)       
  t_p_II_de <- rbind(t_p_II_de, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_II_de_Other_cor"]]


t_r_II_de_df <- as.data.frame(t_r_II_de)%>%
  rownames_to_column(var = "row")

write_xlsx(t_r_II_de_df,"~/Downloads/seq/aab/cor/aab_II_de_r.xlsx")  


t_p_II_de_df <- as.data.frame(t_p_II_de)%>%
  rownames_to_column(var = "row")

write_xlsx(t_p_II_de_df,"~/Downloads/seq/aab/cor/aab_II_de_p.xlsx")  





###pearson
library(Hmisc)
# Get all dataset names that match the pattern
datasets <- ls(pattern = "^pivoted_II_de_.*_cor$")

# Create an empty list to store correlation results
correlation_results <- list()
t_r_II_de_pear <- NULL
t_p_II_de_pear <- NULL
# Loop through each dataset
for (dataset in datasets) {
  # Get the data
  data <- get(dataset)
  
  # Compute Spearman correlation
  correlation <- rcorr(as.matrix(data), type = "pearson")
  r_data <- correlation$r 
  p_data <- correlation$P
  # Store the result in the list, using dataset name as the key
  correlation_results[[dataset]] <- correlation
  
  assign(paste0(dataset, "_r"), r_data)
  assign(paste0(dataset, "_p"), p_data)
  
  t_r_II_de_pear <- rbind(t_r_II_de_pear, r_data)       
  t_p_II_de_pear <- rbind(t_p_II_de_pear, p_data)       
}

# Access correlation results by dataset name
# Example: correlation_results[["pivoted_II_de_Other_cor"]]


t_r_II_de_df_pear <- as.data.frame(t_r_II_de_pear)%>%
  rownames_to_column(var = "row")

write_xlsx(t_r_II_de_df_pear,"~/Downloads/seq/aab/cor/aab_II_de_r_pear.xlsx")  


t_p_II_de_df_pear <- as.data.frame(t_p_II_de_pear)%>%
  rownames_to_column(var = "row")

write_xlsx(t_p_II_de_df_pear,"~/Downloads/seq/aab/cor/aab_II_de_p_pear.xlsx")  





