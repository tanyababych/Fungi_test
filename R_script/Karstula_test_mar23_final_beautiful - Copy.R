library(dplyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(data.table)
library(rdist)
library(adespatial)
library(devtools)
library(factoextra)
library(ggpubr)
library(stringr)
library(tibble)
library(ExpDes)
library(tidyverse)
library(RColorBrewer)

#set working directory
setwd(dir = "C:/Users/martit/Desktop/Karstula_Tania_test/")
getwd()

#load("CWD_alldata.RData")


#load the data

otutab = read.delim("karstula_fungi_otu_table.txt", sep = "\t", row.names = 1)
colSums(otutab)
rowSums(otutab)
min(rowSums(otutab))

#loading taxonomy table
tax = read.delim("karstula_fungi_unite8.2_blastn_best_corr.csv", sep = ",",  header = TRUE, row.names = 1)


mapping = read.delim("metadata_karstula.csv", sep = ",", row.names = 1)

#merge otu and tax to filter
otutax = left_join(rownames_to_column(otutab, "OTU_ID"), rownames_to_column(tax, "OTU_ID"), by = "OTU_ID")

#filtering

otutax_filtered = otutax %>% 
    transform(Genus=str_replace(Genus,"g__",""),
            Phylum=str_replace(Phylum, "p__", ""),
            Kingdom=str_replace(Kingdom, "k__", "")) %>%
  filter(Counts>1 & Kingdom == "Fungi") 


#%>% 
 # mutate(Genus = case_when(Similarity...<90 ~ "Unclassified", 
                           #TRUE ~ Genus))

unique(otutax_filtered$Kingdom)



##Fungal ecology
ecology = read.delim("ecology.csv", sep = ",") #I am introducing ecology table separately
#Add ecology and filter
otutax_filtered$Ecology = ecology$primary_lifestyle[match(otutax_filtered$Genus, ecology$Genus_corrected)]

n_distinct(otutax_filtered$Ecology) #checking how many different ecologies I have in my dataset

unique(otutax_filtered$Ecology) 

otutax_filtered_eco = otutax_filtered %>%
  mutate(Ecology = ifelse(is.na(Ecology), "", Ecology)) %>% 
  mutate(Ecology_new = case_when(
    str_detect(Ecology, 'saprotroph') ~ "Saprotroph",
    str_detect(Ecology, 'plant_pathogen') ~ "Plant pathogen",
    str_detect(Ecology, 'ectomycorrhizal') ~ "EcM",
    str_detect(Ecology, "arbuscular") ~ "AM",
    str_detect(Ecology, 'unspecified') ~ "Unknown",
    str_detect(Ecology, "") ~ "Others", 
    TRUE ~ "Unknown"
  ))

unique(otutax_filtered_eco$Ecology_new)



#Add Genus_Phylum column
otutax_filtered_eco$Genus_phylum = paste0(otutax_filtered_eco$Genus, " (",otutax_filtered_eco$Phylum,")")
head(unique(otutax_filtered_eco$Genus_phylum)) 


otutax_filtered_eco$Phylum_genus = paste0(otutax_filtered_eco$Phylum, " (",otutax_filtered_eco$Genus,")")


#split again into filtered otu and tax (to use in phyloseq)
otutab_filtered = otutax_filtered_eco %>% 
  select(1:49) %>% 
  remove_rownames %>% 
  column_to_rownames(var="OTU_ID")

tax_filtered = otutax_filtered_eco %>% 
  select(OTU_ID, Kingdom:Species, Ecology:Phylum_genus) %>% 
  remove_rownames %>% 
  column_to_rownames(var="OTU_ID")

colnames(tax_filtered)
#phyloseq objects

otu_fin = otu_table(otutab_filtered, taxa_are_rows = TRUE)
tax_fin = tax_table(as.matrix(tax_filtered))
mapping_fin = sample_data(mapping)

finmer = phyloseq(otu_fin, tax_fin, mapping_fin)




#alpha diversity


min_depth = min(colSums(otutab_filtered)) 
t_otutable = as.data.frame((t(otutab_filtered))) #transpose otu table for use in vegan, rows are samples, columns are otus

t_otus_rarefied = as.data.frame(round(rrarefy(t_otutable, min_depth)))

?rrarefy
?specnumber

#Calculating richness, shannon, evenness
richness = as.data.frame(specnumber(t_otus_rarefied, MARGIN = 1)) %>% setnames("richness") #number of OTUs (=molecular species) in each sample
shannon = as.data.frame(diversity(t_otus_rarefied, index = "shannon", MARGIN = 1, base = exp(1))) %>% setnames("shannon") #shannon diversity
evenness = as.data.frame(shannon/log(richness)) %>% setnames("evenness") #pielous evenness 

head(richness) #head means top 10 lines 
head(shannon)
head(evenness)

#combining all these indices into one table 
diversity = print(list(rownames_to_column(mapping, "Sample_ID"), 
                       rownames_to_column(richness, "Sample_ID"), 
                       rownames_to_column(shannon, "Sample_ID"), 
                       rownames_to_column(evenness, "Sample_ID"))) %>%
  reduce(full_join, by = "Sample_ID") %>% na.omit()


max(diversity$richness) 

#statistical summary 

diversity_summary = diversity  %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
            mean_richness= mean(richness),
            sd_richness = sd(richness),
            n_richness = n(), 
            SE_richness = sd(richness)/sqrt(n()),
            mean_shannon= mean(shannon),
            sd_shannon = sd(shannon),
            n_shannon = n(), 
            SE_shannon = sd(shannon)/sqrt(n()),
            mean_evenness= mean(evenness),
            sd_evenness = sd(evenness),
            n_evenness = n(), 
            SE_evenness = sd(evenness)/sqrt(n()))





#Compute the analysis of variance
richness.aov = aov(richness ~ Treatment, data = diversity)
shannon.aov = aov(shannon ~ Treatment, data = diversity)
evenness.aov = aov(evenness ~ Treatment, data = diversity)

# Summary of the analysis
summary(richness.aov)
summary(shannon.aov)
summary(evenness.aov)

#No significance 
crd(diversity$Treatment, diversity$richness, quali = TRUE, mcomp = "tukey", sigF = 0.05)
crd(diversity$Treatment, diversity$shannon, quali = TRUE, mcomp = "tukey", sigF = 0.05)



plot_richness = ggplot(diversity_summary, aes(Treatment, mean_richness)) + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = mean_richness - SE_richness, ymax = mean_richness + SE_richness), 
                width = 0.2) +
  theme_classic() +
  theme(
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 14)) +
  labs( x= "Treatment", y = "Taxa richness (mean)") #+
plot_richness


plot_shannon = ggplot(diversity_summary, aes(Treatment, mean_shannon), fill = "white") + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = mean_shannon - sd_shannon, ymax = mean_shannon + sd_shannon), 
                width = 0.2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 14)) +
  labs( x= "Treatment", y = "Shannon divrsity")#+
plot_shannon

plot_evenness = ggplot(diversity_summary, aes(Treatment, mean_evenness), fill = "white") + 
  geom_col(color = "black", fill = "white") +
  geom_errorbar(aes(ymin = mean_evenness - sd_evenness, ymax = mean_evenness + sd_evenness), 
                width = 0.2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 14)) +
  labs( x= "Treatment", y = "Taxa evenness") #+
plot_evenness

ggarrange(plot_richness, plot_shannon, plot_evenness, nrow = 1)



#ORdination - beta diversity 
finmer

finmer.ord = ordinate(finmer, "NMDS", "bray")
nmds = plot_ordination(finmer, finmer.ord, type="samples", color="Treatment") 
nmds

nmds + geom_polygon(aes(fill=Treatment)) + geom_point(size=5) + ggtitle("samples")


#Ordination with vegan package 

otu_hell = as.data.frame(decostand(t_otutable, "hell")) #hellinger transformation (sqrt of relative abundance), this solves rare taxa issue
head(otu_hell)



###NMDS
mapping2 = rownames_to_column(mapping, "Sample_ID") 
veganNMDS = metaMDS(otu_hell, distance = "bray") #metaMDs function from vegan, bray-curtis distance is used, this is similar to ordinate from phyloseq
veganNMDS #note the final stress, this should be reported on graph or figure caption. If higher than 0.2 it is usually not a good indicator of groups
NMDS1 = veganNMDS$points[,1] #taking coordinates from the vegan object
NMDS2 = veganNMDS$points[,2]

#fixing coordinate names
NMDS1_corr = as.data.frame(NMDS1) #coordinates to data.frame, so that we can merge it
setDT(NMDS1_corr, keep.rownames = TRUE)[] #rownames to first column 
colnames(NMDS1_corr) = c("Sample_ID", "NMDS1") #renaming columns to match mappping (metadata)

NMDS2_corr = as.data.frame(NMDS2)
setDT(NMDS2_corr, keep.rownames = TRUE)[]
colnames(NMDS2_corr) = c("Sample_ID", "NMDS2")

mergeNMDS = full_join(NMDS1_corr, NMDS2_corr, by = "Sample_ID") #combine NMDS coordinates
final_NMDS = full_join(mergeNMDS, mapping2, by = "Sample_ID", na.rm = TRUE) #combine coordinates and mapping by "Sample_ID"

ordiplot(veganNMDS, display = "sites") #plot NMDS basic
#noticed an obvious outlier
plot(veganNMDS, display = "sites", type = "text") 

#define custom colors, we have 8 treatments, so we need 8 colors

color_custom = c("burlywood3", "orange", "darkgreen", "darkred", "#673770", "#D14285",  "aquamarine4","#652926")

#Plot without outlier
NMDS_plot = ggplot(final_NMDS, aes(x=NMDS1, y=NMDS2, colour = Treatment)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title= "Fungal community NMDS") +
  labs( col = "Treatment") +
  theme(legend.text=element_text(size=15),
        legend.position = "right",
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
   scale_color_manual(values = color_custom) 

NMDS_plot

#Relative abundance, using phyloseq

finmer

genus.sum_1st_ngs = tapply(taxa_sums(finmer), tax_table(finmer)[, "Genus"], sum, na.rm=TRUE) #sum counts per each genus to see which are most abundant
phylum.sum_1st_ngs = tapply(taxa_sums(finmer), tax_table(finmer)[, "Phylum"], sum, na.rm=TRUE) #sum counts per each phylum to see which are most abundant

#I take a look how many total genera or phyla I have

#top X genera, if you want just phyla ignore and put Phylum in the for loop 



topgenus = names(sort(genus.sum_1st_ngs, TRUE))[1:30] #I decided to have a look into top 30 genera, and everything else will be "other fungi"
other_genus = names(sort(genus.sum_1st_ngs, TRUE))[31:458] #I have 458 genera total

topphylum = names(sort(phylum.sum_1st_ngs, TRUE))[1:3] #I will have a look into top 3 phyla
other_phylum = names(sort(phylum.sum_1st_ngs, TRUE))[4:19] #the rest of phyla are other fungi

#switch to relative abundance (still phyloseq object)
finmer_relabund = transform_sample_counts(finmer,function(x) {x/sum(x)}) #this is all otus



#prune taxa: trim the dataset to the most abundant taxa (number chosen in previous step)
#note the loss of OTUs, showing really just a portion of the data

finmer_topgenus = prune_taxa((tax_table(finmer_relabund)[, "Genus"] %in% topgenus), finmer_relabund)
finmer_others_genus = prune_taxa((tax_table(finmer_relabund)[, "Genus"] %in% other_genus), finmer_relabund)


#Checking number of OTUs, number of topgenus + othergenus should be same as original
finmer

finmer_topgenus
finmer_others_genus

#use psmelt to get the data from phyloseq objects 

finmer_topgenus_melted = psmelt(finmer_topgenus)
finmer_others_genus_melted = psmelt(finmer_others_genus)

colnames(finmer_topgenus_melted)

#Naming all genus as Other fungi in others
finmer_others_genus_melted$Genus_phylum = "Other fungi"
finmer_others_genus_melted$Phylum_genus = "Other fungi"

#replace unclassified with Other fungi in both top genus melted file and top phylum melted file
finmer_topgenus_melted_corr = transform(finmer_topgenus_melted, Genus_phylum = replace(Genus_phylum, startsWith(Genus_phylum, "unclassified"), "Other fungi"))

#Now merge top genus with other genus
finmer_topgenus_others_combined = rbind(finmer_topgenus_melted_corr, finmer_others_genus_melted)

n_distinct(finmer_topgenus_others_combined$Genus_phylum)

colnames(finmer_topgenus_others_combined)

#tax_level_G = c("Genus") 
#tax_level_P = c("Phylum")
#tax_level_F = c("Family")


#make tables by aggregating data that you need for plots (only the metadata you need etc) 
head(finmer_topgenus_others_combined)

#finmer_genus_phylum_merged


finmer_aggregate_all = aggregate(finmer_topgenus_others_combined$Abundance, 
                                 by=list(Samples=finmer_topgenus_others_combined$Sample, 
                                         OTU_ID = finmer_topgenus_others_combined$OTU,
                                         Genus = finmer_topgenus_others_combined$Genus,
                                         Genus_phylum = finmer_topgenus_others_combined$Genus_phylum,
                                         Phylum= finmer_topgenus_others_combined$Phylum,
                                         Phylum_genus = finmer_topgenus_others_combined$Phylum_genus,
                                         Treatment = finmer_topgenus_others_combined$Treatment, 
                                         Block_ID = finmer_topgenus_others_combined$BlockID,
                                         Ecology = finmer_topgenus_others_combined$Ecology,
                                         Ecology_new = finmer_topgenus_others_combined$Ecology_new), FUN = sum)

Relative_abundance = (prop.table(finmer_aggregate_all$x)*100) #relative abundance fix, making proportion to 100 
sum(Relative_abundance) #check point, should be 100

finmer_genus_fin = cbind(finmer_aggregate_all, Relative_abundance) 
#finmer_phylum_fin = cbind(finmer_aggregate_all_P, Relative_abundance)


#add proper colnames if needed

colnames(finmer_genus_fin) = c("Samples", "OTU_ID","Genus","Genus_phylum", "Phylum", "Phylum_genus" , "Treatment", "Block_ID", "Ecology","Ecology_new","Abundance", "Relative_abundance")



#PLOTTING RELATIVE ABUNDANCE


#put other fungi first for pretty plotting :)
finmer_genus_fin$Genus_phylum = relevel(factor(finmer_genus_fin$Genus_phylum), "Other fungi") 
finmer_genus_fin$Phylum_genus = relevel(factor(finmer_genus_fin$Phylum_genus), "Other fungi") 
finmer_genus_fin$Ecology_new = relevel(factor(finmer_genus_fin$Ecology_new), "Others") 


finmer_genus_fin$Ecology_new = factor(finmer_genus_fin$Ecology_new, levels = c("Others", "Unknown", "EcM","AM", "Plant pathogen", "Saprotroph"))

# some more pretty colors

testcolors = c("#eeeeee",
              "burlywood3", "antiquewhite4",
               "cornsilk1", "chocolate4",
               "lightgoldenrod", "lemonchiffon3",
               "navajowhite4", "sandybrown",
               "tan3", "tan4",
               "rosybrown4", "wheat4",
               "wheat1", "seashell4",
               "tan", "sienna", 
              "#003300", "chartreuse4",
               "#339900", "darkgreen",
               "darkolivegreen", "darkseagreen4",
               "olivedrab", "olivedrab4",
               "mediumseagreen", "springgreen3",
               "springgreen4", "palegreen4",
               "darkseagreen2","paleturquoise4", 
               "paleturquoise4", "lightgreen")

testcolors2 = c("#eeeeee" ,"burlywood3", "#065535", "#99CCFF", "darkseagreen4", 
           "#3333CC", "#660033", "#FF9966",  "#CBD588", "#5F7FC7", 
           "orange","#DA5724", "#508578", "#CD9BCD", "darkolivegreen4", "darkred" ,"#673770",
           "#D14285",  "aquamarine4","#652926",
           "#8569D5", "#5E738F")



phylum = print (ggplot(finmer_genus_fin, aes(x=Treatment, y = Relative_abundance, fill = Phylum)) + 
                  theme_bw() +
                  geom_bar(stat='identity', position = "fill", width = .8) + #bar width can be set here
                  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size  = 14, face = "bold"),
                        axis.text.y = element_text(size = 14, face = "bold"),
                        axis.title.x = element_blank()) +
                  scale_fill_manual(values = testcolors) +
                  labs( x= "Decomposition stage", y = "Relative abundance") +
                  guides(fill = guide_legend(ncol = 1)) +
                  theme(strip.text.x = element_text(size = 14, face = "bold"),
                        legend.text = element_text(size = 18, face = "italic"),
                        legend.position = "right",
                        axis.title.y = element_text(size = 14, face = "bold"),
                        axis.title.x = element_text(size = 14, face = "bold"), 
                        legend.title = element_text(size = 14)))

topphyla = finmer_genus_fin[finmer_genus_fin$Phylum %in% topphylum,]
otherphyla = finmer_genus_fin[finmer_genus_fin$Phylum %in% other_phylum,]
otherphyla$Phylum = "Other fungi"
topphyla_corr = transform(topphyla, Phylum = replace(Phylum, startsWith(Phylum, "unclassified"), "Other fungi"))


phylum_combined = rbind(topphyla_corr, otherphyla) %>% arrange(Phylum)
phylum_combined = phylum_combined[order(phylum_combined$Phylum),]
class(phylum_combined)

phylum_combined$Phylum = relevel(factor(phylum_combined$Phylum), "Other fungi") 


phylum_top = print (ggplot(phylum_combined, aes(x=Treatment, y = Relative_abundance, fill = Phylum)) + 
                      theme_bw() +
                      geom_bar(stat='identity', position = "fill", width = .8) + #bar width can be set here
                      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size  = 14, face = "bold"),
                            axis.text.y = element_text(size = 14, face = "bold"),
                            axis.title.x = element_blank()) +
                      scale_fill_manual(values = c("#eeeeee" ,"burlywood3", "darkgreen", "darkred", "darkseagreen4", 
                                                   "tan4")) +
                      labs( x= "Decomposition stage", y = "Relative abundance") +
                      guides(fill = guide_legend(ncol = 1)) +
                      theme(strip.text.x = element_text(size = 14, face = "bold"),
                            legend.text = element_text(size = 18, face = "italic"),
                            legend.position = "right",
                            axis.title.y = element_text(size = 14, face = "bold"),
                            axis.title.x = element_text(size = 14, face = "bold"), 
                            legend.title = element_text(size = 14))) #+
#facet_grid(.~Depth, scales = "free_x", drop = TRUE))

phylum_top

                
                
genus = 
  print (ggplot(finmer_genus_fin, aes(x=Treatment, y = Relative_abundance, fill = Genus_phylum)) + 
                 theme_bw() +
                 geom_bar(stat='identity', position = "fill", width = .8) + #bar width can be set here
                 theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size  = 14, face = "bold"),
                       axis.text.y = element_text(size = 14, face = "bold"),
                       axis.title.x = element_blank()) +
           scale_fill_manual(values = testcolors2) +
                    labs( x= "Treatment", y = "Relative abundance") +
                 guides(fill = guide_legend(ncol = 1)) +
                 theme(strip.text.x = element_text(size = 14, face = "bold"),
                       legend.text = element_text(size = 18, face = "italic"),
                       legend.position = "right",
                       axis.title.y = element_text(size = 14, face = "bold"),
                       axis.title.x = element_text(size = 14, face = "bold"), 
                       legend.title = element_text(size = 14))) #+
                 #facet_grid(.~Depth, scales = "free_x", drop = TRUE))
head(finmer_genus_fin)

genus


ecology = print (ggplot(finmer_genus_fin, aes(x=Treatment, y = Relative_abundance, fill = Ecology)) + 
                   theme_bw() +
                   geom_bar(stat='identity', position = "fill", width = .8) + #bar width can be set here
                   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size  = 14, face = "bold"),
                         axis.text.y = element_text(size = 14, face = "bold"),
                         axis.title.x = element_blank()) +
                   scale_fill_manual(values = c("#eeeeee" ,"burlywood3", "#065535", "#99CCFF", "darkseagreen4", 
                                                "#3333CC", "#660033", "#FF9966",  "#CBD588", "#5F7FC7", 
                                                "orange","#DA5724", "#508578", "#CD9BCD", "darkolivegreen4", "darkred" ,"#673770",
                                                "#D14285",  "aquamarine4","#652926",
                                                "#8569D5", "#5E738F")) +
                   labs( x= "", y = "Relative abundance") +
                   guides(fill = guide_legend(ncol = 1)) +
                   theme(strip.text.x = element_text(size = 14, face = "bold"),
                         legend.text = element_text(size = 18, face = "italic"),
                         legend.position = "right",
                         axis.title.y = element_text(size = 14, face = "bold"),
                         axis.title.x = element_text(size = 14, face = "bold"), 
                         legend.title = element_text(size = 14)) +
                   facet_grid(.~Treatment, scales = "free_x", drop = TRUE))
ecology


  ecology = print (ggplot(finmer_genus_fin, aes(x=Treatment, y = Relative_abundance, fill = Ecology_new)) + 
                   theme_bw() +
                   geom_bar(stat='identity', position = "fill", width = .8) + #bar width can be set here
                   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size  = 14, face = "bold"),
                         axis.text.y = element_text(size = 14, face = "bold"),
                         axis.title.x = element_blank()) +
                   scale_fill_manual(values = c("#eeeeee" ,"burlywood3", "#065535", "#99CCFF", "darkseagreen4", 
                                                "#3333CC", "#660033", "#FF9966",  "#CBD588", "#5F7FC7", 
                                                "orange","#DA5724", "#508578", "#CD9BCD", "darkolivegreen4", "darkred" ,"#673770",
                                                "#D14285",  "aquamarine4","#652926",
                                                "#8569D5", "#5E738F")) +
                   labs( x= "", y = "Relative abundance") +
                   guides(fill = guide_legend(ncol = 1)) +
                   theme(strip.text.x = element_text(size = 14, face = "bold"),
                         legend.text = element_text(size = 18, face = "italic"),
                         legend.position = "right",
                         axis.title.y = element_text(size = 14, face = "bold"),
                         axis.title.x = element_text(size = 14, face = "bold"), 
                         legend.title = element_text(size = 14)) +
                   facet_grid(.~Treatment, scales = "free_x", drop = TRUE))
  


