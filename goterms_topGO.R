library("topGO")
library("Rgraphviz")

setwd("/proj/uppstore2017185/b2014034_nobackup/Aurora/Migratory_divide/Enrichment_TopGO")

# 1. These steps will lod our data and create our GO annotation
geneID2GO <- readMappings(file = "Vanessa_TopGo_reference_base.tsv")

geneUniverse <- names(geneID2GO)

genesOfInterest.bv <- read.table("input_goterms_fst_outlier_gene_names_025.txt", header=FALSE)
genesOfInterest.bv <- unlist(genesOfInterest.bv)
genesOfInterest.bv <- as.character(genesOfInterest.bv)


geneList.bv <- factor(as.integer(geneUniverse %in% genesOfInterest.bv))
names(geneList.bv) <- geneUniverse

# 2. Now we make the topgo data object and choose the GO hierarchy: BP for Biological Process, MF for Molecular Function or CC for Cellular Component.
myGOdata.bv.bp <- new("topGOdata", description="Candidate genes", ontology="BP", allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize= 5)
myGOdata.bv.mf <- new("topGOdata", description="Candidate genes", ontology="MF", allGenes=geneList.bv,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize= 5)

# 3. Test for significance. We have to define the type of test.
## Note: Selecting ‘algorithm=classic’ means that the GO hierarchy isn’t taken into account, so each GO term is tested independently (over-representation enrichment). 
## The limitation of using this is that all genes annotated to a GO terms will be automatically annotated to its parents as well, therefore a GO term might look enriched just because its children are enriched. 
## Thus, it is important that GO hierarchy is taken into account (conditional enrichment) to avoid redundancy.
## This means using weight01
classic.Fisher.bv.bp <- runTest(myGOdata.bv.bp, algorithm="classic", statistic="fisher")
weight.Fisher.bv.bp <- runTest(myGOdata.bv.bp, algorithm="weight", statistic="fisher")
weight01.Fisher.bv.bp <- runTest(myGOdata.bv.bp, algorithm="weight01", statistic="fisher")
classic.Fisher.bv.mf <- runTest(myGOdata.bv.mf, algorithm="classic", statistic="fisher")
weight.Fisher.bv.mf <- runTest(myGOdata.bv.mf, algorithm="weight", statistic="fisher")
weight01.Fisher.bv.mf <- runTest(myGOdata.bv.mf, algorithm="weight01", statistic="fisher")

# Examine and explore the data: extract pvalues from results using weight algorithm, classic, and see which are in common
pvalFis.mf <- score(weight.Fisher.bv.mf)
pvalFis.bp <- score(weight.Fisher.bv.bp)

pValue.classic.mf <- score(classic.Fisher.bv.mf)
pValue.weight.mf <- score(weight.Fisher.bv.mf)[names(pValue.classic.mf)]

gstat.mf <- termStat(myGOdata.bv.mf, names(pValue.classic.mf))
gSize.mf <- gstat.mf$Annotated / max(gstat.mf$Annotated) * 4

pValue.classic.bp <- score(classic.Fisher.bv.bp)
pValue.weight.bp <- score(weight.Fisher.bv.bp)[names(pValue.classic.bp)]

gstat.bp <- termStat(myGOdata.bv.bp, names(pValue.classic.bp))
gSize.bp <- gstat.bp$Annotated / max(gstat.bp$Annotated) * 4

# 4. Obtain the results -> this is without multiple correction!!!
allRes.bv.bp <- GenTable(myGOdata.bv.bp,
                         classic = classic.Fisher.bv.bp,
                         weight = weight.Fisher.bv.bp,
                         weight01 = weight01.Fisher.bv.bp,
                         orderBy = "classic", ranksOf = "classic", topNodes = 200)

allRes.bv.mf <- GenTable(myGOdata.bv.mf,
                         classic = classic.Fisher.bv.mf,
                         weight = weight.Fisher.bv.mf,
                         weight01 = weight01.Fisher.bv.mf,
                         orderBy = "classic", ranksOf = "classic", topNodes = 200)



# 5. Correct for multiple testing
## performing BH correction on our p values for the weight algorithm
p.adj.bp=round(p.adjust(allRes.bv.bp$weight01,method="BH"),digits = 4)
p.adj.mf=round(p.adjust(allRes.bv.mf$weight01,method="BH"),digits = 4)

p.adj.bp.fdr=round(p.adjust(allRes.bv.bp$weight01,method="fdr"),digits = 4)
p.adj.mf.fdr=round(p.adjust(allRes.bv.mf$weight01,method="fdr"),digits = 4)


## create the file with all the statistics from GO analysis
allRes.bv.bp.corrBH=cbind(allRes.bv.bp,p.adj.bp)
allRes.bv.mf.corrBH=cbind(allRes.bv.mf,p.adj.mf)
allRes.bv.bp.corrBH=allRes.bv.bp.corrBH[order(allRes.bv.bp.corrBH$p.adj.bp),]
allRes.bv.mf.corrBH=allRes.bv.mf.corrBH[order(allRes.bv.mf.corrBH$p.adj.mf),]

allRes.bv.bp.corrFDR=cbind(allRes.bv.bp,p.adj.bp.fdr)
allRes.bv.mf.corrFDR=cbind(allRes.bv.mf,p.adj.mf.fdr)
allRes.bv.bp.corrFDR=allRes.bv.bp.corrFDR[order(allRes.bv.bp.corrFDR$p.adj.bp.fdr),]
allRes.bv.mf.corrFDR=allRes.bv.mf.corrFDR[order(allRes.bv.mf.corrFDR$p.adj.mf.fdr),]



## get list of significant GO before multiple testing correction
results.table.bp.bef= allRes.bv.bp.corrBH[which(allRes.bv.bp.corrBH$weight01<=0.001),]
results.table.mf.bef= allRes.bv.mf.corrBH[which(allRes.bv.mf.corrBH$weight01<=0.001),]

results.table.bp.bef= allRes.bv.bp.corrFDR[which(allRes.bv.bp.corrFDR$weight01<=0.001),]
results.table.mf.bef= allRes.bv.mf.corrFDR[which(allRes.bv.mf.corrFDR$weight01<=0.001),]

## get list of significant GO after multiple testing correction
results.table.bp.corrBH=allRes.bv.bp.corrBH[which(allRes.bv.bp.corrBH$p.adj.bp<=0.05),]
results.table.mf.corrBH=allRes.bv.mf.corrBH[which(allRes.bv.mf.corrBH$p.adj.mf<=0.05),]


results.table.bp.corrFDR=allRes.bv.bp.corrFDR[which(allRes.bv.bp.corrFDR$p.adj.bp<=0.05),]
results.table.mf.corrFDR=allRes.bv.mf.corrFDR[which(allRes.bv.mf.corrFDR$p.adj.mf<=0.05),]


## save the ontolgies sorted by adjusted pvalues
write.table(allRes.bv.bp.corrBH[1:50,],"Fst_outliers_corrBH_50.csv",sep=",",quote=FALSE,row.names=FALSE) # The first 50

write.table(allRes.bv.bp.corrBH,"TopGO_FST_outliers_025_corrBH_bp.csv",sep=",",quote=FALSE,row.names=FALSE) 
write.table(allRes.bv.mf.corrBH,"TopGO_FST_outliers_025_corrBH_mf.csv",sep=",",quote=FALSE,row.names=FALSE)

write.table(allRes.bv.bp.corrFDR,"TopGO_FST_outliers_025_corrFDR_bp.csv",sep=",",quote=FALSE,row.names=FALSE) 
write.table(allRes.bv.mf.corrFDR,"TopGO_FST_outliers_025_corrFDR_mf.csv",sep=",",quote=FALSE,row.names=FALSE)

## Examples of Dasha's way:
write.table(allRes.bv.bp[,c(1,6)], file = "Branch_TopGO_revigo_bp.tsv", sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes.bv.bp, file = "Branch_TopGO_go_bp.tsv", sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE)
write.table(allRes.bv.mf[,c(1,6)], file = "Branch_TopGO_revigo_mf.tsv", sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(allRes.bv.mf, file = "Branch_TopGO_go_mf.tsv", sep = "\t", qmethod = "double", quote = FALSE, row.names = FALSE)


# 6. Make a visualization plot. Change the number of significant nodes to the desired number
pdf(file='topGOPlot_FSToutliers_bp.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(myGOdata.bv.bp, score(weight.Fisher.bv.bp), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()

# Dasha's way
printGraph(myGOdata.bv.bp, classic.Fisher.bv.bp, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


# 7.Get all genes annotated in your significant GOterms
myterms = results.table.bp.corrBH$GO.ID
myterms = results.table.mf.corrBH$GO.ID
mygenes = genesInTerm(myGOdata.bv.bp)

var=c("number", "GOID", "geneID")
for (i in 1:length(myterms))
{
myterm=myterms[i]
mygenesforterm= mygenes[myterm][[1]]
mygenesforterm=paste(mygenesforterm, collapse=',')
var[i]=paste(myterm,"\t",mygenesforterm)
}


gene2GO_significant <- read_table2(var)
colnames(gene2GO) <- c("GOID", "geneID")


write.table(gene2GO_significant,"genetoGOmapping.bp.025.prova.txt", sep="\t", row.names=FALSE)
write.table(gene2GO_significant,"genetoGOmapping.mf.025.txt", sep="\t", quote=F)

# 8. Now we subset the genes of interest (those in the initial list=input_goterms_fst_outlier_gene_names_025.txt)
mygenes = genesInTerm(myGOdata.bv.bp) # load full annotation 
anotation = lapply(mygenes,function(x) x[x %in% genesOfInterest.bv] ) # select the annotations for my genes of interest 
threshold <- 0.05 # this is the threshold of the p-values of the GO terms
significantGO <- subset(allRes.bv.bp.corrBH, p.adj.bp < threshold) # select significant GO terms
gene_interest_GO_significant <- anotation[significantGO$GO.ID] # intersect the genes of interest with significant GOterms

gene_interest_GO_significant_df <- data.frame(lapply(tableGO.genes, "length<-", max(lengths(tableGO.genes)))) # convert it to a dataframe
gene_interest_GO_significant_dft <- t(gene_interest_GO_significant_df) # transpose it and it's ready
write.table(gene_interest_GO_significant_dft,"genetoGOmapping.bp.025.onlygenesofinterest.txt", col.names=FALSE, sep="\t")


# Repeat for molecular function



