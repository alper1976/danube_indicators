
#### load second chunk of libraries####
if (!require("vegan")) {
  install.packages("vegan", dependencies = TRUE)
  library(vegan)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("ape")) {
  install.packages("ape", dependencies = TRUE)
  library(ape)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)
###load PHYLOSEQ####
#BiocManager::install("phyloseq")


library(phyloseq)

if (!require("Hmisc")) {
  install.packages("Hmisc", dependencies = TRUE)
  library(Hmisc)
}
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}
if (!require("Imap")) {
  install.packages("Imap", dependencies = TRUE)
  library(Imap)
}
if (!require("MASS")) {
  install.packages("MASS", dependencies = TRUE)
  library(MASS)
}
if (!require("ade4")) {
  install.packages("ade4", dependencies = TRUE)
  library(ade4)
}
if (!require("clustsig")) {
  install.packages("clustsig", dependencies = TRUE)
  library(clustsig)
}
if (!require("matrixStats")){
  install.packages("matrixStats")
  library(matrixStats)
}
if (!require("gridExtra")){
  install.packages("gridExtra")
  library(gridExtra)
}
if (!require("stats")){
  install.packages("stats")
  library(stats)
}
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(ggrepel)
}
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
}
if (!require("picante")){
  install.packages("picante", dependencies = TRUE)
  library(picante)
}
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("mice")){
  install.packages("mice")
  library(mice)
}

if (!require("R.utils")) {
  install.packages("R.utils", dependencies = TRUE)
  library(R.utils)
}
if (!require("pls")) {
  install.packages("pls", dependencies = TRUE)
  library(pls)
}
if (!require("mdatools")) {
  install.packages("mdatools", dependencies = TRUE)
  library(mdatools)
}

if (!require("EcoSimR")) {
  install.packages("EcoSimR", dependencies = TRUE)
  library(EcoSimR)
}
if (!require("Biostrings")) {
  install.packages("Biostrings", dependencies = TRUE)
  library(Biostrings)
}
if (!require("RgoogleMaps")) {
  install.packages("RgoogleMaps", dependencies = TRUE)
  library(RgoogleMaps)
}
if (!require("dada2")) {
  install.packages("dada2", dependencies = TRUE)
  library(dada2)
}
if (!require("cowplot")) {
  install.packages("cowplot", dependencies = TRUE)
  library(cowplot)
}
if (!require("plsdepot")) {
  install.packages("plsdepot", dependencies = TRUE)
  library(plsdepot)
}
if (!require("readxl")) {
  install.packages("readxl", dependencies = TRUE)
  library(readxl)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
if (!require("ggrepel")) {
  install.packages("ggrepel", dependencies = TRUE)
  library(ggrepel)
}






##################### 3 --- J3S  #Survey ####


#qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual",] 

#color_palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color_palette <- c("firebrick","gray10","dodgerblue","darkolivegreen",
                   "deeppink","steelblue4","tomato","gold1",
                   "red4","cyan4","grey25","navy",
                   "dodgerblue4","palegreen4","indianred1",
                   #v2colors:
                   "gray69","yellow2", "blue4","violetred2",
                   "steelblue1", "darkorange3","deeppink4","chartreuse",
                   "coral1","aquamarine","midnightblue",
                   #AlexEi:
                   "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
                   "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466",
                   "#882255", "#AA4499")




#### ASV table in analysis
rotu <- as.data.frame (t(asv.j3s))

dim(asv.j3s)
#ma <- cbind(ma, "observation"=1:nrow(ma)) #add another column and then rename the rownames

getSequences(asv.j3s)

str(rotu)
variables<-md.j3s


names(variables)[names(variables) == "PicoGreen_DNA-concentration__[ng/µL]"] <- "PicoGreen_DNA_concentration"
names(variables)[names(variables) == "bed_material_D16__[mm]"] <- "bed_material_D16"
names(variables)[names(variables) == "bed_material_D50__[mm]"] <- "bed_material_D50"
names(variables)[names(variables) == "bed_material_D84__[mm]"] <- "bed_material_D84"
names(variables)[names(variables) == "impounded_TRUE-FALSE"] <- "impounded"
names(variables)[names(variables) == "Log+1(AllBac)__[log+1(ME/100mL)]"] <- "Log1_AllBac"
names(variables)[names(variables) == "Log+1(Total_Coliforms)__[log+1(MPN/100ml)]"] <- "Log1_Total_Coliforms"
names(variables)[names(variables) == "Log+1(E.coli)__[log+1(MPN/100ml)]"] <- "Log1_E.coli"
names(variables)[names(variables) == "Log+1(BacR)__[log+1(ME/100mL)]"] <- "Log1_BacR"
names(variables)[names(variables) == "Log+1(BacHum)__[log+1(ME/100mL)]"] <- "Log1_BacHum"
names(variables)[names(variables) == "Log+1(Enterococci)__[log+1(MPN/100mL)]"] <- "Log1_Enterococci"
names(variables)[names(variables) == "Log+1(HF183II___[log+1(ME/100mL)]_[ME/100mL])`"] <- "Log1_HF183II"
names(variables)[names(variables) == "Log+1(Pig2Bac)__[log+1(ME/100mL)]"] <- "Log1_Pig2Bac"
names(variables)[names(variables) == "Ntot__[mg/L]"] <- "Ntot"
names(variables)[names(variables) == "Carbon__NUMERIC__[ng/L]__DETECTION-LIMIT=20"] <- "Carbon"
names(variables)[names(variables) == "Ptot__[mg/L]"] <- "Ptot"
names(variables)[names(variables) == "DOC_mg/L"] <- "DOC"
names(variables)[names(variables) == "Chl-a_Phytoplankton_CHECK_comparability_JDS2vs3"] <- "Chl_a"


###INSERT WATER QUALITY CLASSES FOR ALL DANUBE SITES, SAME VALUE FOR LEFT, RIGHT , MIDDLE####
####Values are taken from Table 7 of JDS3 report, these values reflect the saprobic index from Multi Habitat Sampling (MHS).  
MHS_WQ<-as.factor(c("Good",
          rep("Good",3), rep("Good", 3), rep("Good", 3),
          rep("Good", 6), rep("Moderate", 3), rep("Good", 12), rep("Good", 3), #10
          rep("Good", 3), "Good", rep("Good", 3), rep("Good", 3), rep("Good", 3), "Good", #16
          rep("Good", 3), "Good", rep("Good", 3), rep("Good", 3), rep("Good", 3), rep("Moderate", 3), #22
          "Good", rep("Good", 3), rep("Good", 3), rep("Good", 3), rep("Good", 3),rep("Poor", 3), #28
          rep("Poor", 3),rep("Moderate", 3), rep("Good", 3), rep("Bad", 3), rep("Good", 3), #33
          rep("Moderate", 3), rep("Moderate", 3), rep("Good", 3),rep("Good", 3),rep("Good", 3), #38
          rep("Moderate", 3), rep("Moderate", 3), "Moderate", rep("Poor", 3),rep("Good", 3), #43
          rep("Good", 3), rep("Moderate", 3), rep("Poor", 3), rep("Good", 3), "Good", rep("Good", 3), #49
          rep("Good", 3), "Good", rep("Good", 3), rep("Good", 3), "Good", rep("High", 3), "High", #56
          rep("High", 3), "High", rep("Good", 3), rep("Good", 3), rep("Moderate", 3), #61
          rep("Good", 3), rep("Good", 3), rep("Moderate", 3), rep("Good", 3), rep("Good", 3) #68 (63 and 64 are missing in the entire dataset and for the reference table)
          ))


variables$MHS_WQ<-MHS_WQ
variables$MHS_WQ
variables$MHS_WQ <- factor(variables$MHS_WQ, levels = c("High", "Good", "Moderate", "Poor", "Bad"))
variables[,c("tributary", "MHS_WQ", "impounded")]

####BOXPLOT VARIABLES####
library(ggpubr)

(bxp <- ggviolin(variables, x = "MHS_WQ", y = "Chl_a", color = "MHS_WQ", add="jitter"))
(bxp3 <- ggviolin(variables, x = "MHS_WQ", y = "Ptot", color = "MHS_WQ", add="jitter"))
(bxp4 <- ggviolin(variables, x = "MHS_WQ", y = "Ntot", color = "MHS_WQ", add="jitter"))
(bxp5 <- ggviolin(variables, x = "MHS_WQ", y = "Carbon", color = "MHS_WQ", add="jitter"))
(bxp6 <- ggviolin(variables, x = "MHS_WQ", y = "Log1_Total_Coliforms", color = "MHS_WQ", add="jitter"))
(bxp7 <- ggviolin(variables, x = "MHS_WQ", y = "Log1_Enterococci", color = "MHS_WQ", add="jitter"))
(bxp8 <- ggviolin(variables, x = "MHS_WQ", y = "Log1_AllBac", color = "MHS_WQ", add="jitter"))

(bxp <- ggboxplot(variables, x = "MHS_WQ", y = "Chl_a", color = "MHS_WQ", add="jitter"))
(bxp3 <- ggboxplot(variables, x = "MHS_WQ", y = "Ptot", color = "MHS_WQ", add="jitter"))
(bxp4 <- ggboxplot(variables, x = "MHS_WQ", y = "Ntot", color = "MHS_WQ", add="jitter"))
(bxp5 <- ggboxplot(variables, x = "MHS_WQ", y = "Carbon", color = "MHS_WQ", add="jitter"))
(bxp6 <- ggboxplot(variables, x = "MHS_WQ", y = "Log1_Total_Coliforms", color = "MHS_WQ", add="jitter"))
(bxp7 <- ggboxplot(variables, x = "MHS_WQ", y = "Log1_Enterococci", color = "MHS_WQ", add="jitter"))
(bxp8 <- ggboxplot(variables, x = "MHS_WQ", y = "Log1_AllBac", color = "MHS_WQ", add="jitter"))

(bxp0 <- ggboxplot(variables,  y = "Water_Temperature", add="jitter", color = "#CC79A7")+
    xlab("WT") + ylab("°C"))
(bxp <- ggboxplot(variables,  y = "Chl_a", add="jitter", color = "#CC79A7")+
    xlab("Chl-a") + ylab("µg/L"))
(bxp3 <- ggboxplot(variables,  y = "Conductivity", add="jitter", color = "#E69F00")+
  xlab("Cond") + ylab("μS/cm"))
(bxp4 <- ggboxplot(variables,  y = "Ntot", add="jitter", color = "#56B4E9")+
  xlab("Ntot") + ylab("mg/L"))
(bxp5 <- ggboxplot(variables,  y = "DOC", add="jitter", color = "#009E73")+
  xlab("DOC") + ylab("mg/L"))
(bxp6 <- ggboxplot(variables,  y = "Log1_Total_Coliforms", add="jitter", color = "#F0E442")+
    xlab("Total Coliformes") + ylab("[log+1(ME/100mL)]"))
(bxp7 <- ggboxplot(variables,  y = "Log1_E.coli", add="jitter", color = "#0072B2")+
    xlab("E.coli") + ylab("[log+1(MPN/100mL)]"))
(bxp8 <- ggboxplot(variables,  y = "Log1_AllBac", add="jitter", color = "#D55E00")+
    xlab("All bacterial cells") + ylab("[log+1(ME/100mL)]"))
(bxp9 <- ggboxplot(variables,  y = "Log1_BacHum", add="jitter", color = "#D55E00")+
    xlab("Bacteria from Human") + ylab("[log+1(ME/100mL)]"))

library(ggpubr)

pdf("Boxplot_variables.pdf",width = 6.60, height = 4.31)
figure1 <- ggarrange(bxp0,bxp,bxp3,bxp4,bxp5,bxp6,bxp7, bxp8,bxp9,
                     labels = c("A", "B", "C", "D","E","F","G","H", "I"),
                     ncol = 3, nrow = 3)
figure1
dev.off()

#### Create a phyloseq object to analyse the microbial community

variables$impounded<- as.factor(variables$impounded)
levels(variables$impounded)[levels(variables$impounded)=="0"] <- "no"
levels(variables$impounded)[levels(variables$impounded)=="1"] <- "yes"

otu_norw<-otu_table(rotu, taxa_are_rows =T)
tax <- tax_table(as.matrix(taxa.all))
raw_physeq <- phyloseq(otu_norw, tax, sample_data(variables))
raw_physeq

#subset phyloseq object removing Eukaryota and Chloroplast####
bac_physeq_clean1 <- subset_taxa(raw_physeq, Kingdom!="Eukaryota")
bac_physeq_clean1
bac_physeq_clean2 <- subset_taxa(bac_physeq_clean1, Order!="Chloroplast")
bac_physeq_clean2
bac_physeq_clean3 <- subset_taxa(bac_physeq_clean2,  Kingdom!="Archaea")
bac_physeq_clean3
physeq <- subset_taxa(bac_physeq_clean3,  Family!="Mitochondria")
physeq
physeq<-subset_samples(physeq, tributary != "1")
physeq
physeq<-prune_taxa(taxa_sums(physeq) > 0, physeq)
physeq
#rarecurve(otu_table(physeq), step=50, cex=0.5)


abund = data.frame(t(abundances(physeq)))
meta.fine = data.frame(sample_data(physeq))
taxonomy_table = data.frame(tax_table(physeq))

rarecurve(t(otu_table(physeq)), step=50, cex=0.5)

# rarefy without replacement####
ps.rarefied = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)), rngseed = 1, replace = FALSE, trimOTUs = TRUE, verbose = TRUE) 
ps.rarefied

sample_sums(physeq)
sample_sums(ps.rarefied)

####Relative abundances####
#fine Phyla####
ps3 <- tax_glom(physeq, "Phylum")
ps3.top <- merge_less_than_top(ps3, top=10)
pdf(file = "relative abundance 10 most abundant phyla.pdf",width = 6.60, height = 4.31)
plot_bar(ps3.top, fill = "Phylum", "Sample") +
  geom_bar(stat="identity") +
  facet_grid(~transect_code.cv,scales= "free_x") +
  scale_fill_manual(values = color_palette)  +
  #coord_flip() +
  ylab("Percentage of Sequences") + ylim(0, 1) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack")
dev.off()

ps3 <- tax_glom(physeq, "Class")
ps3.top <- merge_less_than_top(ps3, top=10)
pdf(file = "relative abundance 10 most abundant class.pdf",width = 6.60, height = 4.31)
plot_bar(ps3.top, fill = "Class", "Sample") +
  geom_bar(stat="identity") +
  facet_grid(~transect_code.cv,scales= "free_x") +
  scale_fill_manual(values = color_palette)  +
  #coord_flip() +
  ylab("Percentage of Sequences") + ylim(0, 1) +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack")
dev.off()

ps3 <- tax_glom(physeq, "Genus")
ps3.top <- merge_less_than_top(ps3, top=10)
pdf(file = "relative abundance 10 most abundant genus.pdf",width = 6.60, height = 4.31)
plot_bar(ps3.top, fill = "Genus", "Sample") +
  geom_bar(stat="identity") +
  facet_grid(~transect_code.cv,scales= "free_x") +
  scale_fill_manual(values = color_palette)  +
  #coord_flip() +
  ylab("Percentage of Sequences") + ylim(0, 1) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "stack")
dev.off()


######### Plots and stats #####################
library(phyloseq)
## Phylum####
taxonomy = "Phylum"

classGlom = tax_glom(physeq, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix') 
tax_table = data.frame(mapply(`/`, data.frame(taxon_table), sample_sums(physeq)) * 100)
rownames(tax_table) = tax_matrix[,taxonomy] 
variants<- apply(tax_table,1,var)
colSums(tax_table[, c(1:160)] > 0)
prevalence_table = tax_table
prevalence_table[prevalence_table > 0] = 1
tax_table$prevalence <- rowSums(prevalence_table)
tax_table$variants <- variants
tax_table
tax_table <- tax_table[order(-tax_table$prevalence, -tax_table$variants),]
rowSums(taxon_table)
dim(taxon_table)
tax_table_clean <- as.matrix(tax_table[1:6, 1:160])
rownames(tax_table_clean)
tr<-t(tax_table_clean)
rownames(tr)<-row.names(metadf.bx)
tax_table_clean<-tr
metadata_fine<-phyloseq::sample_data(classGlom)[,c(140,431)]
metadata_fine$impounded<-recode_factor(metadata_fine$impounded, "0" = "NO", 
                                               "1" = "Yes")

t_tax<-merge(tax_table_clean, metadata_fine, by="row.names")

melted_tax_table = reshape2::melt(t_tax)
names(melted_tax_table)[1] <- "Var2"
names(melted_tax_table)[4] <- "Var1"
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))

class_colors <- setNames(color_palette, levels(melted_tax_table$Var1))

pdf("Phyla_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", alpha = .5) +
  guides(fill = guide_legend(title = taxonomy)) +
  coord_flip() +
  theme(axis.text = element_text(size=5),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "right",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey")) +
  xlab("samples") +
  ylab("proportion of reads [%]") +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(limits = rev(levels(melted_tax_table$Var2)))+labs(title = "Phyla for JDS3 Survey")
dev.off()

## RDA

metadata4 = data.frame(sample_data(classGlom))

names(metadata4)[names(metadata4) == "PicoGreen_DNA-concentration__[ng/µL]"] <- "PicoGreen_DNA_concentration"
names(metadata4)[names(metadata4) == "bed_material_D16__[mm]"] <- "bed_material_D16"
names(metadata4)[names(metadata4) == "bed_material_D50__[mm]"] <- "bed_material_D50"
names(metadata4)[names(metadata4) == "bed_material_D84__[mm]"] <- "bed_material_D84"
names(metadata4)[names(metadata4) == "impounded_TRUE-FALSE"] <- "impounded"
names(metadata4)[names(metadata4) == "Log+1(AllBac)__[log+1(ME/100mL)]"] <- "Log1_AllBac"
names(metadata4)[names(metadata4) == "Log+1(Total_Coliforms)__[log+1(MPN/100ml)]"] <- "Log1_Total_Coliforms"
names(metadata4)[names(metadata4) == "Log+1(E.coli)__[log+1(MPN/100ml)]"] <- "Log1_E.coli"
names(metadata4)[names(metadata4) == "Log+1(BacR)__[log+1(ME/100mL)]"] <- "Log1_BacR"
names(metadata4)[names(metadata4) == "Log+1(BacHum)__[log+1(ME/100mL)]"] <- "Log1_BacHum"
names(metadata4)[names(metadata4) == "Log+1(Enterococci)__[log+1(MPN/100ml)]"] <- "Log1_Enterococci"
names(metadata4)[names(metadata4) == "Log+1(HF183II___[log+1(ME/100mL)]_[ME/100mL])`"] <- "Log1_HF183II"
names(metadata4)[names(metadata4) == "Log+1(Pig2Bac)__[log+1(ME/100mL)]"] <- "Log1_Pig2Bac"
names(metadata4)[names(metadata4) == "Ntot__.mg.L."] <- "Ntot"
names(metadata4)[names(metadata4) == "Carbon__.ng.L."] <- "Carbon"
names(metadata4)[names(metadata4) == "Ptot__.mg.L."] <- "Ptot"
names(metadata4)[names(metadata4) == "DOC_mg.L"] <- "DOC"
names(metadata4)[names(metadata4) == "Chl.a_Phytoplankton_CHECK_comparability_JDS2vs3"] <- "Chl_a"
names(metadata4)[names(metadata4) == "TSS_CHECK_comparability_JDS2vs3"] <- "TSS"


metadata4 = metadata4[, c("pH", "Conductivity", "Water_Temperature", "Ntot", "Chl_a", "DOC", "Log1_Total_Coliforms", "Log1_BacHum","Log1_AllBac", "Log1_E.coli", "River_km_Distance_to_mouth")]

#"River_km_Distance_to_mouth",
#### check dataset for NAs content and fill missing values ####

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(metadata4,2,pMiss)
apply(metadata4,1,pMiss)
library(mice)
md.pattern(metadata4)
if (!require("VIM")) {
  install.packages("VIM", dependencies = TRUE)
  library(VIM)
}
aggr_plot <- aggr(metadata4, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(metadata4), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

tempData <- mice(metadata4,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
tempData$imp$Ntot
tempData$imp$Log1_BacHum
tempData$imp$Chl_a
tempData$imp$Conductivity
tempData$imp$Water_Temperature
tempData$imp$DOC

completedData <- complete(tempData,1)
metadata4 <- completedData
#permanova####

library(vegan)
library(MASS)
library(stats)
library(ggplot2)
library(ggrepel)



#standardization of the ASV's abundances by using the hellinger method
raOTU_01 <- decostand(tax_table_clean, method = "hellinger") 



###VARIATION PARTITIONING STEP
set.seed(123)
rda1chem<-rda(raOTU_01~1,metadata4)###model with intercept only
rda2chem<-rda(raOTU_01~.,metadata4)###model with all explanatory variables
step.forward<-ordistep(rda1chem,scope=formula(rda2chem),
                       direction="forward",perm.max=200,pstep=999)

anova(step.forward,by="margin", permutations = 999)
#best model from forward selection
raOTU_01 ~ Chl_a + River_km_Distance_to_mouth + pH + Log1_Total_Coliforms +      DOC + Log1_AllBac + Conductivity
#### Multicollinearity test####
if (!require("olsrr")) {
  install.packages("olsrr", dependencies = TRUE)
  library(olsrr)
}

library(oslrr)

testdata<-metadata4[,c("Chl_a", "River_km_Distance_to_mouth","pH", "Log1_Total_Coliforms","DOC","Log1_AllBac", "Conductivity")] #on the basis of forward selection model

sysvalmodel<- lm(Chl_a~River_km_Distance_to_mouth + pH + Log1_Total_Coliforms + DOC + Log1_AllBac + Conductivity, data = testdata)

summary(sysvalmodel)

ols_coll_diag(sysvalmodel)

round(cor(testdata),2)
library(vegan)
library(grid)
abund_table.adonis <- adonis(raOTU_01 ~ ., data=testdata)
abund_table.adonis

bioenv(wisconsin(raOTU_01)~ Chl_a + pH + Log1_Total_Coliforms +DOC + Log1_AllBac + Conductivity, data=metadata4)

###rda-edcostanda with hellinger transform ###
rda_div_meta <- rda(decostand(tax_table_clean, method="hellinger"), as.matrix(metadata4[,c("Chl_a","pH", "Log1_Total_Coliforms","Log1_AllBac", "Conductivity","River_km_Distance_to_mouth")]))

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors)

mult = 0.2

pdf("RDA_Phyla_JDS3Survey.pdf",width = 6.60, height = 4.31)
rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                                                      y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey", alpha = .8) +
  geom_text_repel(data = data.frame(rda_scores_env),
                  aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                      hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
                  color = "grey", size = 6, alpha = .8) +
  # coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.3, 0.2)) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")
rda_biplot_class
dev.off()
rda_biplot_class


### scatterplot 

class_colors_clean<- class_colors[1:6]
n_color <- nrow(tax_table)-length(class_colors_clean)
class_colors_total<- c(class_colors_clean, rep("black",n_color))
names(class_colors_total)= rownames(tax_table)

pdf("PREV_VAR_Phyla_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(tax_table, aes(x=prevalence, y= variants, color = rownames(tax_table))) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors_total)+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")+
  ylab("variance")
dev.off()


## Class####
taxonomy = "Class"

classGlom = tax_glom(physeq, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix') 
tax_table = data.frame(mapply(`/`, data.frame(taxon_table), sample_sums(physeq)) * 100)
rownames(tax_table) = tax_matrix[,taxonomy] 
variants<- apply(tax_table,1,var)

colSums(tax_table[, c(1:160)] > 0)

prevalence_table = tax_table

prevalence_table[prevalence_table > 0] = 1

tax_table$prevalence <- rowSums(prevalence_table)

tax_table$variants <- variants

tax_table

tax_table <- tax_table[order(-tax_table$prevalence, -tax_table$variants),]

rowSums(taxon_table)
dim(taxon_table)
tax_table_clean <- as.matrix(tax_table[1:9, 1:160])
rownames(tax_table_clean)
tr<-t(tax_table_clean)
rownames(tr)<-row.names(metadata4)
tax_table_clean<-tr
metadata_fine<-phyloseq::sample_data(classGlom)[,c(140,431)]
metadata_fine$impounded<-recode_factor(metadata_fine$impounded, "0" = "NO", 
                                       "1" = "Yes")

t_tax<-merge(tax_table_clean, metadata_fine, by="row.names")

melted_tax_table = reshape2::melt(t_tax)
names(melted_tax_table)[1] <- "Var2"
names(melted_tax_table)[4] <- "Var1"
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))

class_colors <- setNames(color_palette, levels(melted_tax_table$Var1))

pdf("Class_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", alpha = .5) +
  guides(fill = guide_legend(title = taxonomy)) +
  coord_flip() +
  theme(axis.text = element_text(size=5),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "right",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey")) +
  xlab("samples") +
  ylab("proportion of reads [%]") +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(limits = rev(levels(melted_tax_table$Var2)))+labs(title = "Class for JDS3 Survey")
dev.off()

## RDA

#permanova####

#standardization of the ASV's abundances by using the hellinger method
raOTU_01 <- decostand(tax_table_clean, method = "hellinger") 

###VARIATION PARTITIONING STEP
rda1chem<-rda(raOTU_01~1,metadata4)###model with intercept only
rda2chem<-rda(raOTU_01~.,metadata4)###model with all explanatory variables
step.forward<-ordistep(rda1chem,scope=formula(rda2chem),
                       direction="forward",perm.max=200,pstep=999)


#best model from forward selection
#raOTU_01 ~ Chl_a + River_km_Distance_to_mouth + pH + Log1_Total_Coliforms + DOC + Log1_AllBac + Conductivity + Water_Temperature 
#### Multicollinearity test####
if (!require("olsrr")) {
  install.packages("olsrr", dependencies = TRUE)
  library(olsrr)
}
library(oslrr)

testdata<-metadata4[,c("Chl_a", "River_km_Distance_to_mouth" , "pH" , "Log1_Total_Coliforms", "DOC" , "Log1_AllBac",  "Conductivity", "Water_Temperature")] #on the basis of forward selection model

sysvalmodel<- lm(Chl_a~River_km_Distance_to_mouth + pH + Log1_Total_Coliforms + DOC + Log1_AllBac + Conductivity + Water_Temperature , data = testdata)

summary(sysvalmodel)

ols_coll_diag(sysvalmodel)

round(cor(testdata),2)
library(vegan)
library(grid)
abund_table.adonis <- adonis(raOTU_01 ~ ., data=testdata)
abund_table.adonis

###rda-edcostanda with hellinger transform ###
rda_div_meta <- rda(decostand(tax_table_clean, method="hellinger"), as.matrix(metadata4[,c("Chl_a", "River_km_Distance_to_mouth" , "pH" , "Log1_Total_Coliforms", "DOC" , "Log1_AllBac",  "Conductivity", "Water_Temperature")]))

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors)

mult = 0.2

pdf("RDA_Class_JDS3Survey.pdf",width = 6.60, height = 4.31)
rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                                                      y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey", alpha = .8) +
  geom_text_repel(data = data.frame(rda_scores_env),
                  aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                      hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
                  color = "grey", size = 6, alpha = .8) +
  # coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.3, 0.2)) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")
rda_biplot_class
dev.off()
rda_biplot_class



### scatterplot 

class_colors_clean<- class_colors[1:9]
n_color <- nrow(tax_table)-length(class_colors_clean)
class_colors_total<- c(class_colors_clean, rep("black",n_color))
names(class_colors_total)= rownames(tax_table)

pdf("PREV_VAR_Class_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(tax_table, aes(x=prevalence, y= variants, color = rownames(tax_table))) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors_total)+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")+
  ylab("variance")
dev.off()


## Genus####
library(ggplot2)
taxonomy = "Genus"

classGlom = tax_glom(physeq, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix') 
tax_table = data.frame(mapply(`/`, data.frame(taxon_table), sample_sums(physeq)) * 100)
rownames(tax_table) = tax_matrix[,taxonomy] 
variants<- apply(tax_table,1,var)

colSums(tax_table[, c(1:160)] > 0)

prevalence_table = tax_table

prevalence_table[prevalence_table > 0] = 1

tax_table$prevalence <- rowSums(prevalence_table)

tax_table$variants <- variants

tax_table

tax_table <- tax_table[order(-tax_table$prevalence, -tax_table$variants),]

rowSums(taxon_table)
dim(taxon_table)
tax_table_clean <- as.matrix(tax_table[1:15, 1:160])
rownames(tax_table_clean)
tr<-t(tax_table_clean)
rownames(tr)<-row.names(metadata4)
tax_table_clean<-tr
metadata_fine<-phyloseq::sample_data(classGlom)[,c(140,431)]
metadata_fine$impounded<-recode_factor(metadata_fine$impounded, "0" = "NO", 
                                       "1" = "Yes")

t_tax<-merge(tax_table_clean, metadata_fine, by="row.names")

melted_tax_table = reshape2::melt(t_tax)
names(melted_tax_table)[1] <- "Var2"
names(melted_tax_table)[4] <- "Var1"
melted_tax_table <- dplyr::arrange(melted_tax_table, Var1, desc(value))

class_colors <- setNames(color_palette, levels(melted_tax_table$Var1))

pdf("Genus_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", alpha = .5) +
  guides(fill = guide_legend(title = taxonomy)) +
  coord_flip() +
  theme(axis.text = element_text(size=5),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "right",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey")) +
  xlab("samples") +
  ylab("proportion of reads [%]") +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(limits = rev(levels(melted_tax_table$Var2)))+labs(title = "Genera for JDS3 Survey")
dev.off()

## RDA
#permanova####

#standardization of the ASV's abundances by using the hellinger method
raOTU_01 <- decostand(tax_table_clean, method = "hellinger") 

###VARIATION PARTITIONING STEP
rda1chem<-rda(raOTU_01~1,metadata4)###model with intercept only
rda2chem<-rda(raOTU_01~.,metadata4)###model with all explanatory variables
step.forward<-ordistep(rda1chem,scope=formula(rda2chem),
                       direction="forward",perm.max=200,pstep=999)


#best model from forward selection
#raOTU_01 ~ River_km_Distance_to_mouth + Water_Temperature + Chl_a +pH + Log1_E.coli + Log1_AllBac + Conductivity

#### Multicollinearity test####
if (!require("olsrr")) {
  install.packages("olsrr", dependencies = TRUE)
  library(olsrr)
}

library(oslrr)

testdata<-metadata4[,c("River_km_Distance_to_mouth" , "Water_Temperature" , "Chl_a" ,      "pH" , "Log1_E.coli" , "Log1_AllBac" , "Conductivity")]#on the basis of forward selection model

sysvalmodel<- lm(River_km_Distance_to_mouth~Water_Temperature + Chl_a +pH + Log1_E.coli + Log1_AllBac + Conductivity
, data = testdata)

summary(sysvalmodel)

ols_coll_diag(sysvalmodel)

round(cor(testdata),2)

library(vegan)
library(grid)
abund_table.adonis <- adonis(raOTU_01 ~ ., data=testdata)
abund_table.adonis
###rda-edcostanda with hellinger transform ###
rda_div_meta <- rda(decostand(tax_table_clean, method="hellinger"), as.matrix(metadata4[,c("River_km_Distance_to_mouth" , "Water_Temperature" , "Chl_a" ,"pH" , "Log1_E.coli" , "Log1_AllBac" , "Conductivity")]))

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

metadata_class = cbind(metadata4, t(taxon_table))

rda_scores_env = vegan::scores(rda_div_meta, display = "bp")
rda_scores_species = vegan::scores(rda_div_meta, display = "sp")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1


rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = RDA1, y = RDA2, color = names)) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors)

mult = 0.2

pdf("RDA_Genus_JDS3Survey.pdf",width = 6.60, height = 4.31)
rda_biplot_class <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * RDA1,
                                                      y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey", alpha = .8) +
  geom_text_repel(data = data.frame(rda_scores_env),
                  aes(x = mult * RDA1, y = mult * RDA2, label = rownames(rda_scores_env),
                      hjust = 0.5 * (1-sign(RDA1)), vjust = 0.5 * (1-sign(RDA2))),
                  color = "grey", size = 6, alpha = .8) +
  # coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.3, 0.2)) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")
rda_biplot_class
dev.off()
rda_biplot_class


### scatterplot 

class_colors_clean<- class_colors[1:15]
n_color <- nrow(tax_table)-length(class_colors_clean)
class_colors_total<- c(class_colors_clean, rep("black",n_color))
names(class_colors_total)= rownames(tax_table)

pdf("PREV_VAR_Genus_JDS3Survey.pdf",width = 6.60, height = 4.31)
ggplot(tax_table, aes(x=prevalence, y= variants, color = rownames(tax_table))) +
  geom_point(size = 10, alpha = .5) +
  scale_color_manual(values = class_colors_total)+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=10),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey"),
        panel.grid = element_line(color = "grey"),
        legend.position = "None")+
  ylab("variance")
dev.off()

