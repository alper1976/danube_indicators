# Bacterial informational content analyses
# Author: Laurent Fontaine

# Load taxonomy table
asv.j3s.taxonomy.filtered = read.table('JS3_taxonomy_filtered.tsv', header=TRUE)

# Read lines in ASV lists file
processFile = function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    print(line)
  }

  close(con)
}

### Chlorophyll a ###

# Load data
dict = read.table("Danube_chl_a_dictionaries.txt", header = TRUE, sep = '\t')
lists = "Danube_chl_a_predictors.txt" # ASVs only

# Top model r-squared values
dict$global_best_r_squared[order(-dict[,"global_best_r_squared"])]

dict[order(-dict[,"global_best_r_squared"]),]
sum(dict$global_best_r_squared > 0.8)/length(dict$global_best_r_squared)

processFile(lists)

dict_chl_a = read.table("Danube_chl_a_dictionaries.txt", header = TRUE, sep = '\t')
dict_chl_a$global_best_r_squared[order(-dict_chl_a[,"global_best_r_squared"])]
dict_chl_a[dict_chl_a$global_best_r_squared == max(dict_chl_a$global_best_r_squared), ]

lists_chl_a = "Danube_chl_a_predictors.txt" # ASVs only
processFile(lists_chl_a)
chl_a_table = read.table(lists_chl_a, sep=',')[dict_chl_a$global_best_r_squared == max(dict_chl_a$global_best_r_squared), ]

# Get a vector of all ASVs in the table to prepare matrix columns
ASV.vec = c()
for (i in chl_a_table){
	ASV.vec = c(ASV.vec, strsplit(i, split='\t')[[1]])
}
ASV.vec

# Create ASV table and fill with presence/absence
ASV.mat.chl_a = as.data.frame(matrix(rep(0, length(chl_a_table)*length(unique(ASV.vec))), length(chl_a_table), length(unique(ASV.vec))))
colnames(ASV.mat.chl_a) = unique(ASV.vec)

temp.vec = c()
for (i in 1:length(chl_a_table)){
	temp.vec = strsplit(as.character(chl_a_table[i]), split='\t')[[1]]
	for (j in temp.vec){
		ASV.mat.chl_a[i,j] = 1
	}
}

# Order ASV matrix by number of occurences
names(colSums(ASV.mat.chl_a))
ASV.mat.chl_a = ASV.mat.chl_a[,order(colSums(-ASV.mat.chl_a))]

# Fix ASV names with correct index to match taxonomy table
output_names = read.table('ML_output_names.tsv', header = TRUE, sep='\t')
name.vec.chl_a = c()
taxonomy.mat.chl_a = matrix(rep("", ncol(ASV.mat.chl_a) * 6), ncol(ASV.mat.chl_a), 6)

for (i in 1:length(colnames(ASV.mat.chl_a))){
	name.vec.chl_a = c(name.vec.chl_a, output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])])
	taxonomy.mat.chl_a[i,1] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[1]])
	taxonomy.mat.chl_a[i,2] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[2]])
	taxonomy.mat.chl_a[i,3] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[3]])
	taxonomy.mat.chl_a[i,4] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[4]])
	taxonomy.mat.chl_a[i,5] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[5]])
	taxonomy.mat.chl_a[i,6] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat.chl_a)[i])],][[6]])
}
colnames(ASV.mat.chl_a) = name.vec.chl_a
rownames(taxonomy.mat.chl_a) = colnames(ASV.mat.chl_a)
colnames(taxonomy.mat.chl_a) = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

chl_a_ASV_vec = colnames(ASV.mat.chl_a)

# Make a matrix with all ASVs found in kept observations and then use presence absence for each ASV for each observation

### Water quality ###
dict_WQ = read.table("Danube_WQ_dictionaries.txt", header = TRUE, sep = '\t')
dict_WQ$global_best_r_squared[order(-dict_WQ[,"global_best_r_squared"])]
dict_WQ[dict_WQ$global_best_r_squared == 1.0, ]

lists_WQ = "Danube_WQ_predictors.txt" # ASVs only
processFile(lists_WQ)
WQ_table = read.table(lists_WQ, sep=',')[dict_WQ$global_best_r_squared == 1.0, ]

# Get a vector of all ASVs in the table to prepare matrix columns
ASV.vec = c()
for (i in WQ_table){
	ASV.vec = c(ASV.vec, strsplit(i, split='\t')[[1]])
}
ASV.vec

# Create ASV table and fill with presence/absence
ASV.mat = matrix(rep(0, length(WQ_table)*length(unique(ASV.vec))), length(WQ_table), length(unique(ASV.vec)))
colnames(ASV.mat) = unique(ASV.vec)

temp.vec = c()
for (i in 1:length(WQ_table)){
	temp.vec = strsplit(as.character(WQ_table[i]), split='\t')[[1]]
	for (j in temp.vec){
		ASV.mat[i,j] = 1
	}
}

# Order ASV matrix by number of occurences
names(colSums(ASV.mat))
ASV.mat = ASV.mat[,order(colSums(-ASV.mat))]

# Fix ASV names with correct index to match taxonomy table
output_names = read.table('ML_output_names.tsv', header = TRUE, sep='\t')
name.vec = c()
taxonomy.mat = matrix(rep("", ncol(ASV.mat) * 6), ncol(ASV.mat), 6)

for (i in 1:length(colnames(ASV.mat))){
	name.vec = c(name.vec, output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])])
	taxonomy.mat[i,1] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[1]])
	taxonomy.mat[i,2] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[2]])
	taxonomy.mat[i,3] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[3]])
	taxonomy.mat[i,4] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[4]])
	taxonomy.mat[i,5] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[5]])
	taxonomy.mat[i,6] = as.character(asv.j3s.taxonomy.filtered[output_names$Correct_names[which(output_names$ML_output_names %in% colnames(ASV.mat)[i])],][[6]])
}
colnames(ASV.mat) = name.vec
rownames(taxonomy.mat) = colnames(ASV.mat)
colnames(taxonomy.mat) = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

# Keep only ASVs occuring in the 90th percentile or above by occurence
ASV.mat = ASV.mat[,colSums(ASV.mat) > quantile(colSums(ASV.mat), 0.95)]

# Get lowest taxonomic level in taxonomy subset
lowest.taxonomic.level.mat = matrix(rep("", nrow(taxonomy.mat)), nrow(taxonomy.mat), 1)
for (i in 1:nrow(taxonomy.mat)){
	if (length(tail(which( !is.na(taxonomy.mat[i, ]), arr.ind=TRUE), n=1))>0){
		lowest.taxonomic.level.mat[i,1] = taxonomy.mat[i, as.numeric(tail(which( !is.na(taxonomy.mat[i, ]), arr.ind=TRUE), n=1))]
	} else {
		lowest.taxonomic.level.mat[i,1] = NA
	}
}
rownames(lowest.taxonomic.level.mat) = rownames(taxonomy.mat)
colnames(ASV.mat) = lowest.taxonomic.level.mat[1:ncol(ASV.mat), 1]

# Generate co-occurrence matrix with crossproduct
library(igraph)
co_mat <- t(ASV.mat) %*% ASV.mat

# Set diagonal values to 0
diag(co_mat) <- 0
 
# Assign dim names
dimnames(co_mat) <- list(colnames(ASV.mat), colnames(ASV.mat))
 
# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
 
# Assign nodes weight equal to species frequency
g <- set.vertex.attribute(g, "v_weight", value = colSums(ASV.mat))
 
pdf('MHS_WQ_ASV_network.pdf')
plot(g, vertex.size = V(g)$v_weight, edge.width = E(g)$weight)
dev.off()

svg('MHS_WQ_ASV_network.svg')
plot(g, vertex.size = V(g)$v_weight, edge.width = E(g)$weight)
dev.off()

# Co-exclusion matrix
co_mat <- t(ASV.mat) %*% (ASV.mat)
co_mat[co_mat>0]=1
co_mat = 1 - co_mat

# Assign dim names
dimnames(co_mat) <- list(colnames(ASV.mat), colnames(ASV.mat))
 
# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
 
# Assign nodes weight equal to species frequency
g <- set.vertex.attribute(g, "v_weight", value = colSums(ASV.mat))
 
pdf('MHS_WQ_co-exclusion_ASV_network.pdf')
plot(g, vertex.size = V(g)$v_weight, edge.width = E(g)$weight)
dev.off()

svg('MHS_WQ_co-exclusion_ASV_network.svg')
plot(g, vertex.size = V(g)$v_weight, edge.width = E(g)$weight)
dev.off()

# Phylogeny
library(phangorn)
library(msa)

name.vec[order(colSums(-ASV.mat))]
name.vec.chl_a
name.vec[order(colSums(-ASV.mat))] %in% name.vec.chl_a
name.vec[order(colSums(-ASV.mat))][name.vec[order(colSums(-ASV.mat))] %in% name.vec.chl_a]

seq.vec.WQ = colnames(asv.j3s)[name.vec[order(colSums(-ASV.mat))]] # Retrieve DNA sequences of top water quality predictors from original JS3 ASV table
seq.vec.chl_a = colnames(asv.j3s)[name.vec.chl_a]
names(seq.vec.WQ) = name.vec[order(colSums(-ASV.mat))]
names(seq.vec.chl_a) = name.vec.chl_a[order(colSums(-ASV.mat.chl_a))]

lowest.taxonomic.level.mat = matrix(rep("", nrow(taxonomy.mat)), nrow(taxonomy.mat), 1)
for (i in 1:nrow(taxonomy.mat)){
	if (length(tail(which( !is.na(taxonomy.mat[i, ]), arr.ind=TRUE), n=1))>0){
		lowest.taxonomic.level.mat[i,1] = paste0(taxonomy.mat[i, as.numeric(tail(which( !is.na(taxonomy.mat[i, 1:3]), arr.ind=TRUE), n=1))], '_', rownames(taxonomy.mat)[i])
	} else {
		lowest.taxonomic.level.mat[i,1] = rownames(taxonomy.mat)[i]
	}
}
rownames(lowest.taxonomic.level.mat) = rownames(taxonomy.mat)

lowest.taxonomic.level.mat.chl_a = matrix(rep("", nrow(taxonomy.mat.chl_a)), nrow(taxonomy.mat.chl_a), 1)
for (i in 1:nrow(taxonomy.mat.chl_a)){
	if (length(tail(which( !is.na(taxonomy.mat.chl_a[i, ]), arr.ind=TRUE), n=1))>0){
		lowest.taxonomic.level.mat.chl_a[i,1] = paste0(taxonomy.mat.chl_a[i, as.numeric(tail(which( !is.na(taxonomy.mat.chl_a[i, 1:3]), arr.ind=TRUE), n=1))], '_', rownames(taxonomy.mat.chl_a)[i])
	} else {
		lowest.taxonomic.level.mat.chl_a[i,1] = rownames(taxonomy.mat.chl_a)[i]
	}
}
rownames(lowest.taxonomic.level.mat.chl_a) = rownames(taxonomy.mat.chl_a)

subset(lowest.taxonomic.level.mat, rownames(lowest.taxonomic.level.mat) %in% name.vec[order(colSums(-ASV.mat))])
subset(lowest.taxonomic.level.mat.chl_a, rownames(lowest.taxonomic.level.mat.chl_a) %in% name.vec.chl_a)

names(seq.vec.WQ) = subset(lowest.taxonomic.level.mat, rownames(lowest.taxonomic.level.mat) %in% name.vec[order(colSums(-ASV.mat))]) #lowest.taxonomic.level.mat
names(seq.vec.chl_a) = subset(lowest.taxonomic.level.mat.chl_a, rownames(lowest.taxonomic.level.mat.chl_a) %in% name.vec.chl_a) #lowest.taxonomic.level.mat.chl_a

!(names(seq.vec.WQ)) %in% names(seq.vec.WQ)[names(seq.vec.WQ) %in%  names(seq.vec.chl_a)]
non_overlapping_ASVs = c(seq.vec.WQ[!(names(seq.vec.WQ)) %in% names(seq.vec.WQ)[names(seq.vec.WQ) %in%  names(seq.vec.chl_a)]], seq.vec.chl_a)

library(seqinr)
write.fasta(as.character(non_overlapping_ASVs), names=names(non_overlapping_ASVs), 'non_overlapping_ASVs.fasta', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(non_overlapping_ASVs, names=names(non_overlapping_ASVs), file.out='non_overlapping_ASVs.fasta', open = "w", as.string = FALSE)

multiple.alignment = msa(non_overlapping_ASVs, method='ClustalW', type='dna', order='input')
seq.vec.phyDAT = as.phyDat(multiple.alignment, type = "DNA", levels = NULL, nbcol=-1)
dna_dist = dist.ml(seq.vec.phyDAT) # K80, GTR nucelotide substitution models
treeNJ = NJ(dna_dist)
fit = pml(treeNJ, data=seq.vec.phyDAT)#, method = "sankoff")
fitGTR = update(fit, k=4, inv=0.2)
fitGTR = optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE, rearrangement='stochastic', control=pml.control(trace=0))
#bts = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, multicore=FALSE, control=pml.control(trace=0))
#plot(midpoint(fitGTR$tree), bts, p=50, type='p')
NJtrees <- bootstrap.phyDat(seq.vec.phyDAT, FUN=function(x)NJ(dist.logDet(x)), bs=100)
treeNJ <- plotBS(treeNJ, NJtrees, "phylogram", p=100) # Mute node bootstrap support

svg('Phylogeny.svg')
plotBS(treeNJ, NJtrees, "phylogram", p=100, cex=0.2)
dev.off()

library(cluster)
library(dendextend)
library(ggplot2)
library(ggjoy)
library(ggtree)
library(dplyr)

taxonomy.df = as.data.frame(taxonomy.mat)
taxonomy.df$Name = names(seq.vec.WQ)
taxonomy.mat.chl_a.df = as.data.frame(taxonomy.mat.chl_a)
taxonomy.mat.chl_a.df$Name = names(seq.vec.chl_a)
taxonomy.df = rbind(taxonomy.df, taxonomy.mat.chl_a.df)
taxonomy.df = taxonomy.df[!duplicated(taxonomy.df), ]
levels(taxonomy.df$Kingdom) = c(levels(taxonomy.df$Kingdom), 'unassigned')
levels(taxonomy.df$Phylum) = c(levels(taxonomy.df$Phylum), 'unassigned')
levels(taxonomy.df$Class) = c(levels(taxonomy.df$Class), 'unassigned')
levels(taxonomy.df$Order) = c(levels(taxonomy.df$Order), 'unassigned')
levels(taxonomy.df$Family) = c(levels(taxonomy.df$Family), 'unassigned')
levels(taxonomy.df$Genus) = c(levels(taxonomy.df$Genus), 'unassigned')
taxonomy.df[is.na(taxonomy.df)] = 'unassigned'
taxonomy.df$Name = as.factor(as.character(taxonomy.df$Name))
rownames(taxonomy.df) = taxonomy.df$Name
taxonomy.dist = daisy(taxonomy.df, metric=c('gower'))
#plot(as.dendrogram(hclust(taxonomy.dist)))
color.vec = c()
symbol.vec = c()
for (i in 1:nrow(taxonomy.df)){
	if (taxonomy.df$Name[i] %in% names(seq.vec.WQ) && taxonomy.df$Name[i] %in% names(seq.vec.chl_a)){
		color.vec = c(color.vec, 'black')
		symbol.vec = c(symbol.vec, 16)
	} else if (taxonomy.df$Name[i] %in% names(seq.vec.WQ)){
		color.vec = c(color.vec, 'white')
		symbol.vec = c(symbol.vec, 15)
	} else {
		color.vec = c(color.vec, 'black')
		symbol.vec = c(symbol.vec, 17)
	}
}
taxonomy.df$color = color.vec
taxonomy.df$symbol = symbol.vec

svg('Phylogeny.svg')
dend <- taxonomy.dist %>% hclust %>% as.dendrogram %>% set('branches_k_color', k=length(unique(taxonomy.df$Class)), value=c('gray', 'blue', 'green', 'yellow', 'magenta', 'brown', 'cyan', 'purple', 'orange', 'pink', 'aquamarine', 'coral', 'dark green', 'orchid', 'red')) %>% 
set("branches_lwd", 0.5) %>% set("labels_cex", c(0.2,0.2)) %>% #set("labels_colors") %>% 
set("leaves_pch", taxonomy.df$symbol) %>% set("leaves_cex", 0.275) %>% set("leaves_col", taxonomy.df$color)
plot(dend, horiz = TRUE)
dev.off()

svg('Phylogeny_legend.svg')
plot(1)
legend('bottomleft', legend=unique(taxonomy.df$Class), fill=c('gray', 'blue', 'green', 'yellow', 'magenta', 'brown', 'cyan', 'purple', 'orange', 'pink', 'aquamarine', 'coral', 'dark green', 'orchid', 'red'), horiz=FALSE)
dev.off()