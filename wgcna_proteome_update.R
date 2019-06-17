library(WGCNA)

pwd <- getwd()
# setwd("WGCNA")

enableWGCNAThreads(nThreads = 2)
options(stringsAsFactors = FALSE)


#=======================================================================#

# Setting WGCNA networkType/TOMtype, default: unsigned
type <- "unsigned"

# Setting whether there is trait data or not, default no trait data => trait_form=0
trait_form <- 0
# trait_form <- 1

# Setting model cutoff of R2, default：0.9
cex1 <- 0.9

#=======================================================================#


# Input Expression data file(蛋白表达谱数据，行：蛋白ID，列：样本名)
expr <- read.table(file = "expression.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
expr <- log2(expr)

# Output1： Protein Expression matrix
res_expr <- data.frame(proteinid = rownames(expr), expr)
write.table(res_expr, file = "Protein.Expr.xls", sep = "\t", col.names = T, row.names = F, quote = F)

# 过滤，蛋白不做这步
# Filter var > 75% genes
# data_filter <- data[which(vars > quantile(vars, probs = seq(0, 1, 0.25))[4]), ]
# dim(data_filter)
# data_filter_expr <- as.data.frame(t(data_filter))

expr <- as.data.frame(t(expr))

# We first check for proteins and samples with too many missing values(这里只删空值大于1/2的proteins，不删样本)
gsg <- goodSamplesGenes(expr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  expr <- expr[, gsg$goodGenes]
}

# Output2： Filter Protein Expression matrix
res_filtered_expr <- data.frame(proteinid = colnames(expr), t(expr))
write.table(res_filtered_expr, file = "Protein.Filtered.Expr.xls", sep = "\t", col.names = T, row.names = F, quote = F)

# Construct a trait data(这里需要区分是性状数据还是样本分组数据)
if (trait_form){
  # 如有性状数据，则读入性状数据(行：样本名，列：性状名，数目要于表达谱数据一致，名字也要一致)
  datTraits <- read.table(file = "trait.txt", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  Traits <- datTraits
}else{
  # 如没有性状数据，则使用样本分组数据，读入样本分组数据（第一列：样本名，第二列：组名）  
  datTraits <- read.table(file = "trait.txt", sep = "\t", header = F, stringsAsFactors = F)
  names(datTraits) <- c("sample", "group")
  group_list <- factor(datTraits$group, levels = unique(datTraits$group))
  Traits <- model.matrix(~group_list+0)
  Traits <- as.data.frame(Traits)
  colnames(Traits) <- levels(group_list)
  rownames(Traits) <- datTraits$sample
}

# Check！如性状数据和表达谱数据不一致，则match一下使性状数据和表达谱数据保持一致
match_traits <- match(row.names(expr), row.names(Traits))
Traits <- Traits[match_traits,]

# Output3： Trait matrix
res_trait <- data.frame(sampleid = rownames(Traits), Traits)
write.table(res_trait, file = "trait.xls", sep = "\t", col.names = T, row.names = F, quote = F)

#========================================================================#
#============================Mudule Directory============================#

if (dir.exists("Module")) {
  setwd("Module")
}else{
  dir.create("Module")
  setwd("Module")
}

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTree <- hclust(dist(expr), method = "average")
png(filename = "Sample_cluster.png", width = 22, height = 15, units = "cm", res = 300)
par(mar = c(2,5,5,2))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex = 0.9, cex.lab = 1.5,
     cex.axis = 1.2, cex.main = 1.2)
dev.off()

# Choosing the soft-thresholding power: analysis of network topology
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(expr, powerVector = powers, networkType = type, verbose = 5)
sft$powerEstimate

png(filename = "Scale_Free_Topology_Model_Fit.png", width = 22, height = 15, units = "cm", res = 300)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Choosing the WGCNA default soft-thresholding power
power <- sft$powerEstimate

#=======================================================================#

# One-step network construction and module detection
net <- blockwiseModules(expr, power = power, maxBlockSize = 7000,
                        networkType = type, minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "WGCNATOM",
                        verbose = 3)

table(net$colors)

mergedColors <- labels2colors(net$colors)

png(filename = "PlotModule.png",width = 18, height = 12, units = "cm", res = 300)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# Output4： Protein ID and expression with module information 
pro2module <- data.frame(proteinid = colnames(expr), module = moduleColors, t(expr))
write.table(pro2module, file = "Module.Protein.xls", sep = "\t", col.names = T, row.names = F, quote = F)
# for(selet_module in unique(moduleColors)){
#   expr_matrix <- t(expr[selet_module == moduleColors])
#   expr_module <- data.frame(proteinid = rownames(expr_matrix), expr_matrix)
#   write.table(expr_module, file = paste0(selet_module, ".Protein.Expr.xls"), sep = "\t", col.names = T, row.names = F, quote = F)
# }

#=======================================================================#

# Quantifying module–trait associations
# calculate module Eigengenes in protein samples(ME矩阵)
MEs0 <- moduleEigengenes(expr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# calculate the corrlation between ME matrix and trait matrix and cor-test
nGenes <- ncol(expr)
nSamples <- nrow(expr)
moduleTraitCor <- cor(MEs, Traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Display correlation values within a heatmap plot
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
png(filename = "Module-trait.relationships.png",width =25, height = 20, units = "cm", res = 300)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=======================================================================#

# Calculate topological overlap anew
dissTOM = 1-TOMsimilarityFromExpr(expr, power = power, corType = "pearson")

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Cprotein the plot function
png(filename = "GeneNetwork_heatmap_plot.png" ,width =25, height = 20, units = "cm", res = 300)
geneTree <- net$dendrograms
# sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all proteins")
dev.off()

setwd(pwd)


#=======================================================================#
#======================Mudule_Enrichment Directory======================#

if (dir.exists("Module_Enrichment")) {
  setwd("Module_Enrichment")
}else{
  dir.create("Module_Enrichment")
  setwd("Module_Enrichment")
}

plyr::d_ply(pro2module, "module", function(x){
  folder <- as.character(unique(x["module"]))
  dir.create(folder)
  setwd(folder)
  
  write.table(x["proteinid"], file = "diff_list.txt", col.names = F, row.names = F, quote = F)
  
  setwd("../")
})

setwd(pwd)

#=======================================================================#
#======================Module_Eigengene Directory=======================#

# Visualizing the network of eigengenes

# if (trait_form == 1){
#   # Combined trait with Module eigengene
#   pwd <- getwd()
#   if (dir.exists("Module_eigengene")) {
#     setwd("Module_eigengene")
#   }else{
#     dir.create("Module_eigengene")
#     setwd("Module_eigengene")
#   }
#   
#   for(select_trait in colnames(Traits)){
#     # Isolate certain trait from the clinical traits
#     select_trait_df <- as.data.frame(Traits[, select_trait])
#     names(select_trait_df) <- select_trait
#     MEs2trait <- orderMEs(cbind(MEs, select_trait_df))
#     
#     
#     png(filename = paste0(select_trait, "_Eigengene.dendrogram.heatmap.png"), width =25, height = 25, units = "cm", res = 300)
#     par(cex = 0.9)
#     plotEigengeneNetworks(MEs2trait, 
#                           setLabels = "", 
#                           marDendro = c(0,8,1,2),
#                           marHeatmap = c(8,8,1,2),
#                           # plotDendrograms = FALSE,
#                           xLabelsAngle = 90)
#     dev.off()
#   }
#   
#   # Just Module eigengene without trait
#   png(filename = "Eigengene.dendrogram.heatmap.png", width =25, height = 25, units = "cm", res = 300)
#   par(cex = 0.9)
#   plotEigengeneNetworks(MEs, 
#                         setLabels = "", 
#                         marDendro = c(0,8,1,2),
#                         marHeatmap = c(8,8,1,2),
#                         # plotDendrograms = FALSE,
#                         xLabelsAngle = 90)
#   dev.off()
#   
# }else{
  
  # Just Module eigengene without trait
  png(filename = "Eigengene.dendrogram.heatmap.png", width =25, height = 25, units = "cm", res = 300)
  par(cex = 0.9)
  plotEigengeneNetworks(MEs, 
                        setLabels = "", 
                        marDendro = c(0,8,1,2),
                        marHeatmap = c(8,8,1,2),
                        # plotDendrograms = FALSE,
                        xLabelsAngle = 90)
  dev.off()
# }

# Output5： Module Eigengenes in each sample
res_ME <- data.frame(sampleid = rownames(expr), MEs)
write.table(res_ME, file = "Module.Eigenproteins.xls", sep = "\t", col.names = T, row.names = F, quote = F)

setwd(pwd)


#=======================================================================#
#====================Protein_relationship Directory=====================#

# Intramodular analysis: identifying genes with high GS and MM

# names (colors) of the modules
modNames <- substring(names(MEs), 3)

# Calculate the correlation between expression matrix and module Eigengenes matrix
geneModuleMembership <- as.data.frame(cor(expr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# names(geneModuleMembership) <- paste("MM", modNames, sep="")
# names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitSignificance <- as.data.frame(cor(expr, Traits, use = "p"))
geneTraitPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))


for(select_trait in colnames(Traits)){
  pwd <- getwd()
  if (dir.exists(paste0("Protein_relationship/",select_trait))) {
    setwd(paste0("Protein_relationship/",select_trait))
  }else{
    dir.create(paste0("Protein_relationship/",select_trait), recursive = T)
    setwd(paste0("Protein_relationship/",select_trait))
  }
  trait_column <- match(select_trait, colnames(Traits))
  
  for(selet_module in modNames){
    module_column <- match(selet_module, modNames)
    
    # extract specified model genes
    moduleGenes = moduleColors == selet_module
    
    # Display the genes highy related with model and trait
    png(filename = paste0(selet_module, ".moduel.", "GS_vs_MM.png"), width =20, height = 15, units = "cm", res = 300)
    
    if (selet_module == "lightcyan"){
      col <- "cyan"
    }else{
      col <- selet_module
    }
    
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitSignificance[moduleGenes, trait_column]),
                       xlab = paste("Module Membership in", selet_module, "module"),
                       ylab = paste("Protein significance for", select_trait),
                       main = paste("Module membership vs. protein significance\n"),
                       cex.main = 1, cex.lab = 1, cex.axis = 1, col = col)
    dev.off()
  }
  setwd(pwd)
}


# Output6： Proteins with GS and MM value
write.table(data.frame(proteinid = rownames(geneModuleMembership), geneModuleMembership), 
            file = "Protein_relationship/Protein_ModuleMembership.xls", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(proteinid = rownames(geneTraitSignificance), geneTraitSignificance), 
            file = "Protein_relationship/Protein_GeneSignificance.xls", sep = "\t", col.names = T, row.names = F, quote = F)


#=======================================================================#
#=====================Protein_ModulePlot Directory======================#

# Displaying module heatmap plot

# selet_module <- "turquoise"
for(selet_module in modNames){
  ME <- MEs[, paste("ME", selet_module, sep="")]
  plotMat_matrix <- t(scale(expr[, moduleColors==selet_module]))
  
  pwd <- getwd()
  if (dir.exists("Protein_ModulePlot")) {
    setwd("Protein_ModulePlot")
  }else{
    dir.create("Protein_ModulePlot")
    setwd("Protein_ModulePlot")
  }
  
  png(filename = paste0(selet_module, ".module.heatmap_barplot.png"), width =20, height = 20, units = "cm", res = 300)
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
  plotMat(plotMat_matrix,
          nrgcols = 30, 
          rlabels = F,
          rcols= selet_module, 
          clabels = colnames(plotMat_matrix), 
          main = "",
          cex.main = 1)
  
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=selet_module, 
          main="", cex.main=2,
          ylab="Eigenprotein expression",
          xlab=paste0(selet_module, " module"))
  
  dev.off()
  
  setwd(pwd)
}


#=======================================================================#
#====================Module_visualization Directory=====================#

# Exporting to Cytoscape
TOM <- TOMsimilarityFromExpr(expr, power = power, corType = "pearson")
probes <- names(expr)

pwd <- getwd()
if (dir.exists("Module_visualization")) {
  setwd("Module_visualization")
}else{
  dir.create("Module_visualization")
  setwd("Module_visualization")
}

for(selet_module in modNames){
  inModule <- is.finite(match(moduleColors, selet_module))
  modProbes <- probes[inModule]
  
  # Select the corresponding Topological Overlap
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(selet_module, collapse="-"), ".txt", sep=""),
                                 # nodeFile = paste("CytoscapeInput-nodes-", paste(selet_module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 # altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}
setwd(pwd)
