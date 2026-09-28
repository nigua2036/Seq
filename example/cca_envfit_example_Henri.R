##CCA Analysis example Henri #####


#Here is the list of all the packages that I've used.
#I'm sorry I can't remember exactly which of these are needed, because I always just load them all!


library("vegan")
library("devtools")
# library("ggvegan")
library("BiodiversityR")
library("BiocManager")
library("dada2")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("car")
library("corrplot")
library("MVN")
library("factoextra")
library("readxl")
library("ggpubr")
library("phyloseq")
library("ggrepel")


install.packages("ggvegan")

#### IMPORTANT!!!
#The dataframes containing the microbial population data, and the other variables need to have the exact same sample names as row names!

#####CCA

# "bottle_vegan_pop_for_NMDS" is a dataframe containing the microbial population data, Rows = samples, Columns = relative abundance of each OTU
bottle_ca <- cca(bottle_vegan_pop_for_NMDS) 



#"bottle_chem_con2_14032023" is a dataframe containing the data for chemical and physical variables. Rows = samples, Columns = data
#The envfit function only puts out vectors for variables for which the P value is less than 0.05
envfit_test_ca_14032023 <- envfit(bottle_ca, bottle_chem_con2_14032023, permu = 999, na.rm = TRUE) 



#This part is about extracting the vectors from the CCA data for use in ggplot graph. ggplot can't otherwise use it (or it is beyond my skills). The rest of this portion is about adjusting the lengts of the vectors so that they look better in the graph.
test14032023_arrows <- as.data.frame(envfit_test_ca_14032023$vectors$arrows)
test14032023_r <- scores(envfit_test_ca_14032023$vectors$r)
test14032023_r <- sqrt(test14032023_r)
test14032023_arrows$CA1 <- test14032023_arrows$CA1 / test14032023_r * 1
test14032023_arrows$CA2 <- test14032023_arrows$CA2 / test14032023_r * 1
test14032023_arrows$varNames <- rownames(test14032023_arrows)




#Trying to extract the bottle_ca scores for species and sites so that I can use the ggplot function because it has better functionality 
bottle_ca_scores <- scores(bottle_ca)
bottle_ca_scores_sites <- as.data.frame(bottle_ca_scores$sites)



#Here I am putting in some data to the data frame so that I can group them in ggplot.
#You can also change the sample names here if they are so large they overlap in the graph. 
bottle_ca_scores_sites$Source_type <- c("Agricultural", "Agricultural", "Agricultural","Agricultural", "Agricultural", "Agricultural","Agricultural", "Agricultural", "Agricultural", "Non-agricultural", "Non-agricultural", "Non-agricultural", "Non-agricultural", "Non-agricultural", "Non-agricultural", "Agricultural", "Agricultural", "Agricultural" )
bottle_ca_scores_sites$Source_type <- factor(bottle_ca_scores_sites$Source_type, levels = c("Agricultural", "Non-agricultural"))

bottle_ca_scores_sites$Sample_source <- c("UAG","UAG","UAG","JAG","JAG","JAG","KAG","KAG","KAG","KFR","KFR","KFR","KAF","KAF","KAF","RAG","RAG","RAG")
bottle_ca_scores_sites$Sample_source <- factor(bottle_ca_scores_sites$Sample_source, levels = c("UAG", "JAG", "RAG", "KAG", "KAF", "KFR"))



#This is the ggplot script
bottle_ca_scores_sites %>%
  ggplot() +
  geom_point(aes(x = CA1, y = CA2, color = Sample_source), size = 5) +
  geom_text(aes(x = CA1, y = CA2, label = Sample_source, vjust = -0.7), check_overlap = TRUE) +
  labs(x = "CA1",
       y = "CA2",
       shape = "Source_type") +
  theme_classic() +
  geom_segment(data = test14032023_arrows,
               mapping = aes(x = 0, y = 0, xend = CA1, yend = CA2), color="blue") +
  geom_label_repel(data = test14032023_arrows,
                   mapping = aes(x = CA1, y = CA2, label = varNames), color="blue") +
  #geom_point(data = bottle_cent, aes(x = NMDS1, y = NMDS2, size = 6)) +
  stat_ellipse(data = bottle_ca_scores_sites, alpha = 0.7, aes(x = CA1, y = CA2, group = Source_type, linetype = Source_type)) +
  labs(x = "CA1",
       y = "CA2",
       linetype = "Soil land use",
       color = "Site") +
  theme_classic() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.90, 0.35),
        legend.background = element_rect(fill = "white",
                                         color = "black"),
        legend.margin = margin(t = 6, r = 2, b = 3, l = 3))







##from the chat
##CCA MOB_all
data.cca <- cca(dat_all ~ Moisture. + CH4ox + CO2prod + org_matter + NH4_cm3 + pH , IDall)
plot(data.cca)

ps<-c(rep(1,9), rep(2,2), rep(5,9), rep(6,9), rep(21,9), rep(24,2), rep(23,9), rep(25,9))
cl<-c(rep(1,9), rep(1,2), rep(1,9), rep(1,9), rep(17,9), rep(17,2), rep(17,9), rep(17,9))

legend(x="topleft", legend=c("DNA-Autumn", "DNA-Winter", "DNA-Spring", "DNA-Summer", "RNA-Autumn", "RNA-Winter", "RNA-Spring", "RNA-Summer"),
       pch=c(1,2,5,6,21,24,23,25), col="black", pt.bg=c(1,1,1,1,17,17,17,17), bty="n", title="Wet", cex=0.9)

plot(vare.dryL3, type ="none") 


NMDS example: 
  #dry vare.dryL1 <- metaMDS(dryL1, trace =... x Henri Siljanen
  

#dry
vare.dryL1 <- metaMDS(dryL1, trace = FALSE)

plot(vare.dryL1, type ="none")

ps_dryL1<-c(rep(1,3), rep(5,3), rep(6,3), rep(21,3), rep(23,3), rep(25,2))
cl_dryL1<-c(rep(1,3), rep(1,3), rep(1,3), rep(17,3), rep(17,3), rep(17,2))

points(vare.dryL1, "sites", pch=ps_dryL1, col="black", bg=cl_dryL1)
legend(x="bottomleft", legend=c("DNA-Autumn", "DNA-Spring", "DNA-Summer", "RNA-Autumn", "RNA-Spring", "RNA-Summer"),
       pch=c(1,5,6,21,23,25), col="black", pt.bg=c(1,1,1,17,17,17), bty="n", title="dry_L1", cex=0.9)


ef_dryL1 <- envfit(vare.dryL1, IDdryL1, permu = 999)
ef_dryL1




###4) What are the abundances of the different AOA clades? (in R)####

####4.1) Importation and parsing of the annotation table in R#####
tax <- read.table("~/Downloads/seq/example/5a-uniques_nochim_match_uchimed_tax_assignments.txt",
                  header = FALSE, sep = "\t")

tax$V3 <- NULL
tax$V4 <- NULL

##Split the taxonomic annotation for each level
library("stringr")
tax <- cbind(tax, str_split_fixed(tax$V2, ";", 11))
tax$V2 <- NULL


##Replace empty cells by NA
tax2 <- as.data.frame(apply(tax, 2, function(x) gsub("^$|^ $", NA, x)))


##Remove columns containing only NA
col_to_remove <- c()

for (col in 1:ncol(tax2)) {
  x <- sum(is.na(tax2[,col]))/nrow(tax2)
  if (x == 1) {
    col_to_remove <- c(col_to_remove, col)
  }
}

if (length(col_to_remove) != 0) {
  tax3 <- tax2[,-col_to_remove]
} else {
  tax3 <- tax2
}

#set column names
names(tax3)[1] <- "ASV"
names(tax3)[-1] <- paste0("l", 1:(ncol(tax3)-1), "_tax")

#Set taxonomic annotations as character variables
for (col in 2:ncol(tax3)) {
  tax3[,col] <- as.character(tax3[,col])
}



#Fill all NAs (by copying the previous taxonomic annotation of each ASV
for (col in 1:ncol(tax3)) {
  for (row in 1:nrow(tax3)) {
    if (is.na(tax3[row,col])) {
      if (!grepl("OTU", tax3[row,col-1]) & !grepl("unassigned", tax3[row,col-1])) {
        tax3[row,col] <- paste0(tax3[row,col-1], "_unassigned")
      } else {
        tax3[row,col] <- tax3[row,col-1]
      }
    }
  }
}


####4.2) Aggregate the ASV table by taxonomic annotations#####
####Merge the ASV and the annotations tables
#asv.tax <- merge(asv3.rel, tax3, by.x = "header", by.y = "ASV")

#check the dimension of the merged table

dim(asv.tax)

#export the table
#write.table(asv.tax, "8c-asv_table2_rel_tax.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#example
asv.tax <- read.table("~/Downloads/seq/example/8c-asv_table2_rel_tax.txt", 
                      header = TRUE, quote = "", sep = "\t")

#Aggregate by the full AOA annotation
l11.tax <- aggregate(asv.tax[,2:19], by=list(annotation=asv.tax$l11_tax), FUN=sum)


l11.tax[1:5,1:10]

#Aggregate by the AOA clade annotation (level 2)
l2.tax <- aggregate(asv.tax[,2:19], by=list(annotation=asv.tax$l2_tax), FUN=sum)
l2.tax[,1:10]


####4.3) Aggregate the ASV table by taxonomic annotations####
# First, define for each replicate which sample it belongs to.


replicates_list <- c("peat02", "peat05", "peat08", "peat17", "peat30",
                     "peat31", "peat32", "peat33", "peat38", "peat39",
                     "peat64", "peat65", "peat66", "peat63", "peat67",
                     "peat69", "peat71", "peat77")

replicates_groups <- c("Kev_BS", "Kev_BS", "Kev_BS", "Kev_VS", "Taz_PP_VS",
                       "Taz_PP_VS", "Taz_PP_BS", "Taz_PP_BS", "Taz_PB_BS", "Taz_PB_BS",
                       "Sei_BS", "Sei_BS", "Sei_BS", "Sei_VS", "Tay_BS",
                       "Tay_BS", "Tay_BS", "Tay_VS")

##transpose the asv.tax
n <- l2.tax$annotation
l2.tax.t <- as.data.frame(t(l2.tax[,-1]))
colnames(l2.tax.t) <- n
l2.tax.t$replicate <- rownames(l2.tax.t)
rownames(l2.tax.t) <- NULL


#adding replicate group
l2.tax.t$sample <- rep(NA, nrow(l2.tax.t))
for (row in 1:nrow(l2.tax.t)) {
  l2.tax.t$sample[row] <- replicates_groups[grep(l2.tax.t$replicate[row], replicates_list)]
}


###
l2.tax.t.mean <- aggregate(l2.tax.t[,1:(ncol(l2.tax.t)-2)], by=list(sample=l2.tax.t$sample), FUN=mean)


l2.tax.t.sd <- aggregate(l2.tax.t[,1:(ncol(l2.tax.t)-2)], by=list(sample=l2.tax.t$sample), FUN=sd)

l2.tax.t.mean_sd <- l2.tax.t%>%
  group_by(sample) %>%
  summarize(mean_D = mean(`NS-Delta`),
            sd_D = sd(`NS-Delta`),
            mean_G = mean(`NS-Gamma`),
            sd_G = sd(`NS-Gamma`),
            mean_Z = mean(`NS-Zeta`),
            sd_Z = sd(`NS-Zeta`),
            N=n(),
            .groups = "drop")


write.table(l2.tax.t.mean, "~/Downloads/seq/example/10a-l2_table_rel_mean.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(l2.tax.t.sd, "~/Downloads/seq/example/10b-l2_table_rel_sd.txt", quote = FALSE, sep = "\t", row.names = FALSE)



###5 bubble plot####
# 5.1) Top-50 ASVs bubble plot
# Parse the label for the bubble plot
asv_num <- sub(";size=[0-9]*;", "", asv.tax$header)
asv_size <- sub(";", "", sub("sq[0-9]*;", "", asv.tax$header))
asv_tax <- sub("_unassigned", "", asv.tax$l11_tax)
asv.tax$label <- paste0(asv_tax, " ", asv_num, " (", asv_size, ")")

asv.tax$label[1:9]


# Keep only the 50-top ASVs
asv.tax$average <- rowMeans(asv.tax[,2:19])
asv.tax.sorted <- asv.tax[order(-asv.tax$average),]

asv.tax.sorted$average[1:20]


# Aggregate the smaller taxonomic bins together
asv.tax.sorted$selection <- rep("discarded", nrow(asv.tax.sorted))
asv.tax.sorted$selection[1:50] <- "retained"
asv.tax.sorted$label[asv.tax.sorted$selection == "discarded"] <- "Other"
asv.tax.sorted$average <- NULL
asv.tax.sorted$selection <- NULL
asv.tax.sorted.top50 <- aggregate(asv.tax.sorted[,2:19], by=list(label=asv.tax.sorted$label), FUN=sum)



# Which fraction of the relative abundance won't be displayed in the buble plot?
mean(as.numeric(asv.tax.sorted.top50[51,-1]))

# Transpose the asv.tax.sorted.top50 dataframe --> into wide? format
n <- asv.tax.sorted.top50$label
asv.tax.sorted.top50.t <- as.data.frame(t(asv.tax.sorted.top50[,-1]))
colnames(asv.tax.sorted.top50.t) <- n
asv.tax.sorted.top50.t$replicate <- rownames(asv.tax.sorted.top50.t)
rownames(asv.tax.sorted.top50.t) <- NULL


# Add the replicates groups to this dataframe
asv.tax.sorted.top50.t$sample <- rep(NA, nrow(asv.tax.sorted.top50.t))
for (row in 1:nrow(asv.tax.sorted.top50.t)) {
  asv.tax.sorted.top50.t$sample[row] <- replicates_groups[grep(asv.tax.sorted.top50.t$replicate[row], replicates_list)]
}

# Compute the mean and for each sample
asv.tax.sorted.top50.t.mean <- aggregate(asv.tax.sorted.top50.t[,1:(ncol(asv.tax.sorted.top50.t)-2)],
                                         by=list(sample=asv.tax.sorted.top50.t$sample),
                                         FUN=mean)

?melt
library(data.table)
# Melt the dataframe and parse it
# Melt the dataframe
molten <- melt(asv.tax.sorted.top50.t.mean, id.vars = "sample")

# Remove null values
molten2 <- molten[molten$value > 0,]

# Create a categorical variable to colour the bubbles
molten2$category <- sub(" .+", "", molten2$variable)
library("ggplot2")

# Re-order the factor variables of the ASV labels
molten2$variable <- factor(molten2$variable, levels = rev(levels(molten2$variable)))
# Compute the bubble_plot
library("ggplot2")

bubble_plot <- ggplot(molten2,aes(sample,variable)) +
  geom_point(aes(size=value, fill=molten2$category),shape=21,color="black") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  ylab("AOA ASVs") +
  xlab("Samples") +
  scale_fill_brewer(palette="Paired", name="AOA Taxonomic\nclade") +
  #scale_fill_discrete(name="Taxonomic\nclade") +
  #scale_fill_manual(values= c("maroon2", "pink", "#000000"), name="Taxonomic\nclade") +
  scale_size(name = "Relative\nabundance")

bubble_plot

svg("plots/asv_bubble_plot.svg", width = 7, height = 8)
bubble_plot
dev.off()



# 5.2) Reference-OTUs bubble plots
# For this plot we will aggregate level-5 taxonomic annotation. We will end up with 3 categories only
l5.tax <- aggregate(asv.tax[,2:19], by=list(annotation=asv.tax$l5_tax), FUN=sum)

# Transpose the dataframe
n <- l5.tax$annotation
l5.tax.t <- as.data.frame(t(l5.tax[,-1]))
colnames(l5.tax.t) <- n
l5.tax.t$replicate <- rownames(l5.tax.t)
rownames(l5.tax.t) <- NULL

# Add the replicates groups to this dataframe
l5.tax.t$sample <- rep(NA, nrow(l5.tax.t))
for (row in 1:nrow(l5.tax.t)) {
  l5.tax.t$sample[row] <- replicates_groups[grep(l5.tax.t$replicate[row], replicates_list)]
}
# Compute the mean and for each sample
l5.tax.t.mean <- aggregate(l5.tax.t[,1:(ncol(l5.tax.t)-2)],
                           by=list(sample=l5.tax.t$sample),
                           FUN=mean)

# Melt the dataframe and parse it
# Melt the dataframe
molten3 <- melt(l5.tax.t.mean, id.vars = "sample")

# Remove null values
molten4 <- molten3[molten3$value > 0,]


# Create a categorical variable to colour the bubble

# Re-order the factor variables of the ASV labels
molten4$variable <- factor(molten4$variable, levels = rev(levels(molten4$variable)))

# Compute the bubble_plot
bubble_plot2 <- ggplot(molten4,aes(sample,variable)) +
  geom_point(aes(size=value, fill=variable),shape=21,color="black") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  ylab("AOA ASVs") +
  xlab("Samples") +
  scale_fill_brewer(palette="Paired", name="AOA Taxonomic\nclade") +
  #scale_fill_discrete(name="Taxonomic\nclade") +
  #scale_fill_manual(values= c("maroon2", "pink", "#000000"), name="Taxonomic\nclade") +
  scale_size(name = "Relative\nabundance")

bubble_plot2

svg("plots/l5_bubble_plot.svg", width = 7, height = 3)
bubble_plot2
dev.off()


###6 NMDS analysis#####
# Transpose the ASV tax table
n <- asv.tax$header
asv.tax.t <- as.data.frame(t(asv.tax[,2:19]))
colnames(asv.tax.t) <- n
asv.tax.t$replicate <- rownames(asv.tax.t)
rownames(asv.tax.t) <- NULL
view(asv.tax.t)

# Add the replicates groups to this dataframe
asv.tax.t$sample <- rep(NA, nrow(asv.tax.t))
for (row in 1:nrow(asv.tax.t)) {
  asv.tax.t$sample[row] <- replicates_groups[grep(asv.tax.t$replicate[row], replicates_list)]
}

# Extract categorial variables from samples names (i.e. the site and the type of vegetation)
asv.tax.t$site <- sub("_[A-Z][A-Z]$", "", asv.tax.t$sample)
asv.tax.t$vegetation <- sub(".*_", "", asv.tax.t$sample)

# Define a color for each ecosystem and horizon
asv.tax.t$site_col <- rep(NA, nrow(asv.tax.t))
asv.tax.t$site_col[asv.tax.t$site == "Kev"] <- "#e41a1c"
asv.tax.t$site_col[asv.tax.t$site == "Taz_PP"] <- "#377eb8"
asv.tax.t$site_col[asv.tax.t$site == "Taz_PB"] <- "#4daf4a"
asv.tax.t$site_col[asv.tax.t$site == "Sei"] <- "#984ea3"
asv.tax.t$site_col[asv.tax.t$site == "Tay"] <- "#ff7f00"

# Define shape for each horizon
asv.tax.t$veg_shp <- rep(NA, nrow(asv.tax.t))
asv.tax.t$veg_shp[asv.tax.t$vegetation == "BS"] <- 1
asv.tax.t$veg_shp[asv.tax.t$vegetation == "VS"] <- 2

# Permanova test: are samples significantly different from each other 
# when grouped by site or by vegetation ?
library("vegan")
str(asv.tax.t)
asv.bray <- as.data.frame(as.matrix(vegdist(asv.tax.t[,1:(ncol(asv.tax.t)-6)], diag = TRUE, upper = TRUE)))

adonis2(asv.bray ~ asv.tax.t$site, asv.bray, permutations = 3000) # p = 0.01466 *
adonis2(asv.bray ~ asv.tax.t$vegetation, asv.bray, permutations = 3000) # p = 0.5938
adonis2(asv.bray ~ asv.tax.t$vegetation, asv.bray, permutations = 3000, 
        strata = asv.tax.t$site) # p = 0.4858

# Run the NMDS
asv.nmds <- metaMDS(asv.tax.t[,1:(ncol(asv.tax.t)-6)], k=2)
asv.nmds



# Stress plot
stressplot(asv.nmds)

# svg("plots/asv_nmds_stressplot.svg", width = 4, height = 4)
stressplot(asv.nmds)
dev.off()
summary(asv.nmds)

# Plot the NMDS
# svg("plots/asv_nmds.svg", width = 5, height = 5)
ordiplot(asv.nmds,type="n")
points(asv.nmds, display = 'sites', col = asv.tax.t$site_col, pch = asv.tax.t$veg_shp, cex = 1.5, lwd = 2)
groupz <- unique(asv.tax.t$site)
colz <- unique(asv.tax.t$site_col)

for(i in seq(groupz)) { 
  ordiellipse(asv.nmds, asv.tax.t$site, kind="se", conf=0.95, label=F, font=1, 
              cex=0.5, col=colz[i], show.groups=groupz[i])
}

legend('topright',title = "Ecosystem:", col=colz, legend=unique(asv.tax.t$site), pch = 16, cex = 0.7, box.col = "white", title.adj = FALSE)

legend('topleft', title = "Vegetation:", legend=unique(asv.tax.t$vegetation), 
       pch = sort(unique(asv.tax.t$veg_shp)), cex = 0.7, box.col = "white", , title.adj = FALSE)
title("ASV NMDS")
dev.off()



ggplot(metadata_nmds,
       aes(x=axis1, y=axis2, color=disease_stat, fill=disease_stat)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.2, show.legend=F) +
  geom_point(show.legend=FALSE) +
  geom_richtext(data=my_legend,
                aes(x=x, y=y, label=label, color=color),
                hjust=0, lineheight=0.75,
                inherit.aes = FALSE, show.legend=FALSE,
                fill = NA, label.color = NA
  ) +
  coord_cartesian(xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8)) +
  labs(title="<span style='color:#999999'>**Healthy individuals**</span> have a different microbiota<br>
       from those with <span style='color:#0000FF'>**diarrhea**</span> and those with diarrhea<br>
       who are <span style='color: #ff0000'><strong>positive for *C. difficile*<strong></span>",
       x="NMDS Axis 1",
       y="NMDS Axis 2",
       caption = "All pairwise comparisions were significant using ADONIS at 0.05 using \nBenjimani-Hochberg correction for multiple comparision") +
  scale_color_manual(name=NULL,
                     breaks=c("NonDiarrhealControl", "DiarrhealControl", "Case"),
                     values=c("gray", "blue", "red"),
                     labels=c("Healthy", "Diarrhea",
                              "*C. difficile* positive")
  )+
  scale_fill_manual(name=NULL,
                    breaks=c("NonDiarrhealControl", "DiarrhealControl", "Case"),
                    values=c("lightgray", "dodgerblue", "pink"),
                    labels=c("Healthy", "Diarrhea",
                             "*C. difficile* positive")
  )+
  theme_classic() +
  theme(
    axis.title.y = element_text(hjust=1),
    axis.title.x = element_text(hjust=0),
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(1.05, 0.95),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.title.position = "plot",
    plot.title = element_markdown(margin=margin(b=1, unit="lines")),
    plot.margin = margin(l=0.5, r=0.5, t=0.5, b=0.5, unit="lines"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot"
  )



?vegan


##vegan packages example####
data(dune)
ord <- decorana(dune)
plot(ord)



ord1 <- decorana(ttaxon_polled_I_replicate_imputed[,2:(ncol(ttaxon_polled_I_replicate_imputed)-4)])

summary(ord1)


#MDS
ord_MDS <- metaMDS(dune, trace = FALSE)
ord_MDS
plot(ord_MDS)

plot(ord_MDS, type = "n")
points(ord_MDS, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord_MDS, display = "species", cex=0.7, col="blue")


plot(ord_MDS, type = "n") |>
  points("sites", cex = 0.8, pch=21, col="red", bg="yellow") |>
  text("species", cex=0.7, col="blue")

data(dune.env)
attach(dune.env)

plot(ord_MDS, disp="sites", type="n")
ordihull(ord_MDS, Management, col=1:4, lwd=3)
ordiellipse(ord_MDS, Management, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord_MDS, Management, col=1:4, draw="polygon")
ordispider(ord_MDS, Management, col=1:4, label = TRUE)
points(ord_MDS, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)



###fitting environmental varibles envfit#####
ord.fit <- envfit(ord_MDS ~ A1 + Management, data=dune.env, perm=999)

cladeI_genus_nmds

clade_env <- ttaxon_polled_I_replicate_imputed%>%
  select("sample", "treatment", "soil", "Cu", "waterlevel")

ord_fit_I <- envfit(cladeI_genus_nmds ~ soil + Cu, data = clade_env, perm = 999)

ordisurf(cladeI_genus_nmds, Cu, add = TRUE)
plot(cladeI_genus_nmds, dis = "site")
plot(ord_fit_I)



#only depression soils
detI

detI_nmds <- metaMDS(detI, k = 2)
stressplot(detI_nmds)

clade_envde <- clade_env%>%
  filter(soil == "de")

ord_detI <- envfit(detI_nmds ~ waterlevel + Cu, data = clade_envde, perm = 999)

ordisurf(detI_nmds, Cu, add = TRUE)
plot(detI_nmds, dis = "site")
plot(ord_detI)



plot(ord_MDS, dis="site")
plot(ord.fit)



#Function ordisurf directly adds a fitted surface onto ordination, but it 
# returns the result of the fitted thinplate spline gam
ordisurf(ord_MDS, A1, add=TRUE)




asv.nmds <- metaMDS(asv.tax.t[,1:(ncol(asv.tax.t)-6)], k=2)

ord.fit_alex <- envfit(asv.nmds ~ site + vegetation, data=asv.tax.t, perm=999)

plot(asv.nmds, dis="site")
plot(ord.fit_alex)



###envfit####
data("varechem")

vare.mds <- metaMDS(varespec, trace = FALSE)
ef <- envfit(vare.mds, varechem, permu = 999)
ef


plot(vare.mds, display  = "sites")
plot(ef, p.max = 0.1)


###Constrained ordination cca####
ord_cca <- cca(dune ~ A1 + Management, data=dune.env)
ord_cca
plot(ord_cca)

anova(ord_cca)

ord_cca_all <- cca(dune ~ ., data=dune.env)
plot(ord_cca_all)


anova(ord_cca_all)

##cladeI

ord_cca_cI <- cca(tI_matrix~ soil + Cu, data = clade_env)

plot(ord_cca_cI)
anova(ord_cca_cI)

#analyse terms separately
anova(ord_cca_cI, by="term", permutations=199)

#by marginal effects (“Type IIIeffects”):
anova(ord_cca_cI, by="mar", permutations=199)

#analyse significance of each axis
anova(ord_cca_cI, by="axis", permutations=199)




##alex cca
asv.tax.t_1<- as.data.frame(as.matrix(asv.tax.t[,1:(ncol(asv.tax.t)-6)], diag = TRUE, upper = TRUE))

adonis2(asv.bray ~ asv.tax.t$site, asv.bray, permutations = 3000) # p = 0.01466 *
adonis2(asv.bray ~ asv.tax.t$vegetation, asv.bray, permutations = 3000) # p = 0.5938
adonis2(asv.bray ~ asv.tax.t$vegetation, asv.bray, permutations = 3000, 
        strata = asv.tax.t$site) # p = 0.4858

ord_cca_alex <- cca(asv.tax.t_1 ~ vegetation , data=asv.tax.t)
plot(ord_cca_alex)


#analyse terms separately
anova(ord_cca, by="term", permutations=199)

#by marginal effects (“Type IIIeffects”):
anova(ord_cca, by="mar", permutations=199)

#analyse significance of each axis
anova(ord_cca, by="axis", permutations=199)


#conditioned or partial ordination
ord_cca_cond <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)

ord_cca_cond


anova(ord_cca_cond, by="term", permutations=499)






