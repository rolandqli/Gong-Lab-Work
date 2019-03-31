options(java.parameters = "-Xmx8000m")
library(edgeR)
library(biomaRt)
library(ggplot2)
library(xlsx)
library(dplyr)

# Import counts table.
table <- read.table("combined_counts_table.txt")
table <- select(table, starts_with("C"), starts_with("half"), starts_with("X1h"), 
                starts_with("X2h"), starts_with("X3"),starts_with("X6"), 
                starts_with("X12"), starts_with("X24"), starts_with("X4"))
table <- table[-1:-4,]
rownames(table) <- unlist(lapply(strsplit(rownames(table), split = ".", fixed = T), function(x)x[1]))
View(table)

# Normalized counts per million. 
normalizedcounts <- cpm(table)
View(normalizedcounts)
write.xlsx2(normalizedcounts, "Normalized_counts.xlsx")

# 10-15 counts matches to CPM of 1.
plot(table[,1], normalizedcounts[,1], ylim = c(0,5), xlim = c(0,20))
abline(v = 10, h = 1)
abline(v = 15, h = 1.75)

# Filtering data.
filtering <- rowSums(normalizedcounts > 1) >= 3
filtered <- filter(table, filtering)
rownames(filtered) <- rownames(as.matrix(which(filtering)))
View(filtered)

# Digital Gene Expression List that holds the counts for each sample and holds the group labels,
# library sizes, dispersion estimates and normalization factors.
dgeobj <- DGEList(filtered)
dgeobj <- calcNormFactors(dgeobj)
dgeobj <- estimateCommonDisp(dgeobj)
dgeobj <- estimateGLMTrendedDisp(dgeobj)
dgeobj <- estimateTagwiseDisp(dgeobj)

# Separating into timepoints.
group <- unlist(lapply(strsplit(colnames(filtered), split = "_"), function(x)x[1]))
group <- factor(group, levels = c("Control","halfhr", "X1hr", "X2hr",
                                  "X3hr","X6hr","X12hr","X24hr","X48hr"))
# Checking similarity between each sample.
pdf("MDSPlot.pdf")
plotMDS(dgeobj, labels = group, cex = .75, xlim= c(-2,4))
dev.off()
                       
# Visualising variation.
pdf("BCVPlot.pdf")
plotBCV(dgeobj)
dev.off()

# Model
design <- model.matrix(~0 + group)
voomobj <- voom(dgeobj, design, plot = T)
fit <- lmFit(voomobj, design)

# Import genes and annotation database.
mart <- useMart("ensembl")
mouse <- useDataset("mmusculus_gene_ensembl",mart=mart)
symbols <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"), 
                 filter = "ensembl_gene_id", values = rownames(filtered), mart = mouse)

# 3hr vs. Control
contrast <- makeContrasts(groupX3hr - groupControl, levels = design)
threehrfit <- contrasts.fit(fit, contrast)
threehr <- eBayes(threehrfit)
threehr <- topTable(threehr, sort.by = "P", n = Inf)
ord <- match(rownames(threehr), symbols[,"ensembl_gene_id"])
threehr <- mutate(threehr, ID = rownames(!!threehr), Symbol = !!symbols[ord, "mgi_symbol"], Description = !!symbols[ord, "description"])
matching <- match(threehr$ID, rownames(filtered))
addons <- select(filtered, starts_with("Control"), starts_with("X3hr"))
threehr <- cbind(threehr, addons[matching,])
controlMean <- rowMeans(select(filtered,starts_with("Control")))
threeMean <- rowMeans(select(filtered, starts_with("X3hr")))
MeanDifference <- threeMean - controlMean 
threehr <- mutate(threehr, controlMeans = controlMean[matching], threeMeans = threeMean[matching])
threehr <- mutate(threehr, MeanDifference = threeMeans - controlMeans)
threehr <- select(threehr, -t, -P.Value, -B)
threehr <- select(threehr, ID, Symbol, everything())
View(threehr)
write.xlsx2(threehr, "3hr vs. Control.xlsx")
pdf("3hr.pdf")
ggplot(threehr, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) + geom_point() + scale_color_manual(breaks = c("FALSE", "TRUE"), values=c("black", "red")) + labs(color = "Significant", x = "Average logCPM", title = "3hr vs. Control")
dev.off()

# 6hr vs. Control
contrast <- makeContrasts(groupX6hr - groupControl, levels = design)
sixhrfit <- contrasts.fit(fit, contrast)
sixhr <- eBayes(sixhrfit)
sixhr <- topTable(sixhr, sort.by = "P", n = Inf)
ord <- match(rownames(sixhr), symbols[,"ensembl_gene_id"])
sixhr <- mutate(sixhr, ID = rownames(!!sixhr), Symbol = !!symbols[ord, "mgi_symbol"], Description = !!symbols[ord, "description"])
matching <- match(sixhr$ID, rownames(filtered))
addons <- select(filtered, starts_with("Control"), starts_with("X6hr"))
sixhr <- cbind(sixhr, addons[matching,])
controlMean <- rowMeans(select(filtered,starts_with("Control")))
sixMean <- rowMeans(select(filtered, starts_with("X6hr")))
MeanDifference <- sixMean - controlMean 
sixhr <- mutate(sixhr, controlMeans = controlMean[matching], sixMeans = sixMean[matching])
sixhr <- mutate(sixhr, MeanDifference = sixMeans - controlMeans)
sixhr <- select(sixhr, -t, -P.Value, -B)
sixhr <- select(sixhr, ID, Symbol, everything())
View(sixhr)
write.xlsx2(sixhr, "6hr vs. Control.xlsx")
pdf("6hr.pdf")
ggplot(sixhr, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) + geom_point() + scale_color_manual(breaks = c("FALSE", "TRUE"), values=c("black", "red")) + labs(color = "Significant", x = "Average logCPM", title = "6hr vs. Control")
dev.off()

# 12hr vs. Control
contrast <- makeContrasts(groupX12hr - groupControl, levels = design)
twelvehrfit <- contrasts.fit(fit, contrast)
twelvehr <- eBayes(twelvehrfit)
twelvehr <- topTable(twelvehr, sort.by = "P", n = Inf)
ord <- match(rownames(twelvehr), symbols[,"ensembl_gene_id"])
twelvehr <- mutate(twelvehr, ID = rownames(!!twelvehr), Symbol = !!symbols[ord, "mgi_symbol"], Description = !!symbols[ord, "description"])
matching <- match(twelvehr$ID, rownames(filtered))
addons <- select(filtered, starts_with("Control"), starts_with("X12hr"))
twelvehr <- cbind(twelvehr, addons[matching,])
controlMean <- rowMeans(select(filtered,starts_with("Control")))
twelveMean <- rowMeans(select(filtered, starts_with("X12hr")))
MeanDifference <- twelveMean - controlMean 
twelvehr <- mutate(twelvehr, controlMeans = controlMean[matching], twelveMeans = twelveMean[matching])
twelvehr <- mutate(twelvehr, MeanDifference = twelveMeans - controlMeans)
twelvehr <- select(twelvehr, -t, -P.Value, -B)
twelvehr <- select(twelvehr, ID, Symbol, everything())
View(twelvehr)
write.xlsx2(twelvehr, "12hr vs. Control.xlsx")
pdf("12hr.pdf")
ggplot(twelvehr, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) + geom_point() + scale_color_manual(breaks = c("FALSE", "TRUE"), values=c("black", "red")) + labs(color = "Significant", x = "Average logCPM", title = "12hr vs. Control")
dev.off()

# 24hr vs. Control
contrast <- makeContrasts(groupX24hr - groupControl, levels = design)
tfhrfit <- contrasts.fit(fit, contrast)
tfhr <- eBayes(tfhrfit)
tfhr <- topTable(tfhr, sort.by = "P", n = Inf)
ord <- match(rownames(tfhr), symbols[,"ensembl_gene_id"])
tfhr <- mutate(tfhr, ID = rownames(!!tfhr), Symbol = !!symbols[ord, "mgi_symbol"], Description = !!symbols[ord, "description"])
matching <- match(tfhr$ID, rownames(filtered))
addons <- select(filtered, starts_with("Control"), starts_with("X24hr"))
tfhr <- cbind(tfhr, addons[matching,])
controlMean <- rowMeans(select(filtered,starts_with("Control")))
tfMean <- rowMeans(select(filtered, starts_with("X24hr")))
MeanDifference <- tfMean - controlMean 
tfhr <- mutate(tfhr, controlMeans = controlMean[matching], tfMeans = tfMean[matching])
tfhr <- mutate(tfhr, MeanDifference = tfMeans - controlMeans)
tfhr <- select(tfhr, -t, -P.Value, -B)
tfhr <- select(tfhr, ID, Symbol, everything())
View(tfhr)
write.xlsx2(tfhr, "24hr vs. Control.xlsx")
pdf("24hr.pdf")
ggplot(tfhr, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) + geom_point() + scale_color_manual(breaks = c("FALSE", "TRUE"), values=c("black", "red")) + labs(color = "Significant", x = "Average logCPM", title = "24hr vs. Control")
dev.off()

# 48hr vs. Control
contrast <- makeContrasts(groupX48hr - groupControl, levels = design)
fehrfit <- contrasts.fit(fit, contrast)
fehr <- eBayes(fehrfit)
fehr <- topTable(fehr, sort.by = "P", n = Inf)
ord <- match(rownames(fehr), symbols[,"ensembl_gene_id"])
fehr <- mutate(fehr, ID = rownames(!!fehr), Symbol = !!symbols[ord, "mgi_symbol"], Description = !!symbols[ord, "description"])
matching <- match(fehr$ID, rownames(filtered))
addons <- select(filtered, starts_with("Control"), starts_with("X48hr"))
fehr <- cbind(fehr, addons[matching,])
controlMean <- rowMeans(select(filtered,starts_with("Control")))
feMean <- rowMeans(select(filtered, starts_with("X48hr")))
MeanDifference <- feMean - controlMean 
fehr <- mutate(fehr, controlMeans = controlMean[matching], feMeans = feMean[matching])
fehr <- mutate(fehr, MeanDifference = feMeans - controlMeans)
fehr <- select(fehr, -t, -P.Value, -B)
fehr <- select(fehr, ID, Symbol, everything())
View(fehr)
write.xlsx2(fehr, "48hr vs. Control.xlsx")
pdf("48hr.pdf")
ggplot(fehr, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) + geom_point() + scale_color_manual(breaks = c("FALSE", "TRUE"), values=c("black", "red")) + labs(color = "Significant", x = "Average logCPM", title = "48hr vs. Control")
dev.off()
