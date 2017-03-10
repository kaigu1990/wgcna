# markdownphotos
library(WGCNA)

#process raw data
data <- read.table(file = "WGCNA_test.csv", sep = ",", stringsAsFactors = F,row.names = 1)
colnames(data) <- data[1,]
data <- data[-1,]
rowname <- rownames(data)
data <- as.data.frame(lapply(data[,1:24], as.numeric))
rownames(data) <- rowname

#delete too much missing data
gsg <- goodSamplesGenes(data, verbose = 3)
gsg$allOK
data <- data[gsg$goodSamples, gsg$goodGenes]

#imputed data analysis
library(mice)
data_m <- mice(data, m=10, method = "pmm", maxit = 5, seed = 1)
data_m <- complete(data_m)

#Calculating the average of each sample
names(data_m) <- unlist(lapply(names(data_m), function(x){
  tmp <- strsplit(x, '\\.\\d+$')[[1]][1]
}))

data_matrix <- as.data.frame(lapply(unique(names(data_m)), function(x){
  sample <- rowMeans(data_m[,which(names(data_m)==x)])})
)
names(data_matrix) <- unique(names(data_m))
data_matrix <- t(data_matrix)

#2. 选择合适的阀值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
# Now choose the "signed" and you can choose "unsigned"
sft = pickSoftThreshold(data_matrix, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#choose the power
sft

#3. 一步法网络构建：One-step network construction and module detection
net = blockwiseModules(data_matrix, power = 24, maxBlockSize = 5000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

#4. 绘画结果展示
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#5.结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

#anlysis without trait
MEs0 <- moduleEigengenes(data, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
#build the sample of diag
datTraits <- diag(8)
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- rownames(data)
colnames(datTraits) <- rownames(data)
nSamples = nrow(data)
#connect sample and MEs
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3.5, 8.5, 3, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#from this example, I choose the magenta module
which.module="magenta";

#calculate eigengenes of each module in samples == MEs
datME <- moduleEigengenes(data, moduleColors)$eigengenes
signif(cor(datME, use = "p"), 2)

#the sign of the correlation between the module eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

# Diagnostics: heatmap plots of module expression
# Diagnostics: displaying module heatmap and the eigengene
par(mfrow=c(2,1), mar=c(1, 2, 4, 1))
ME=datME[, paste("ME",which.module, sep="")]
plotMat(t(scale(data[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

# We have a module membership value for each gene in each module
datKME=signedKME(data, datME, outputColumnName="MM.")

#get kME of each gene
a <- as.matrix(datKME)
data_kong <- data.frame(a_name=c(),b_name=c())
for (i in 1:4849) {
  a <- datKME[i,]
  a_name <- rownames(a)
  a <- unlist(as.list(a))
  b <- paste("MM",moduleColors[i],sep = ".")
  b_name <- names(a[b])
  
  e <- data.frame(a_name, b_name, a[b])
  data_kong <- rbind(data_kong,e)
}
data_module <- data_kong[order(data_kong[,2]),]

write.csv(data_kong, file = "protein2module.csv", quote = F)
write.csv(t(data),file = "expression.csv",quote = F)


#1. 可视化全部基因网络
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(data, power = 12);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
#sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#6. 导出网络到Cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(data_matrix_mv, power = 6);
# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
modules = c("turquoise");
# Select module probes选择模块探测
probes = colnames(data_matrix_mv)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

edge <- read.table(file = "AS-green-FPKM-One-step-CytoscapeInput-edges-red.txt", sep = "\t", stringsAsFactors = F)

