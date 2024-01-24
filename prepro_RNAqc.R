##################### Quality control y bias removal ###########################
# PREPROCESAMIENTO

# Packages ---------------------------------------------------------------------

# BiocManager::install("SummarizedExperiment")
# BiocManager::install("biomaRt")
# install.packages ("vroom")
# BiocManager::install("NOISeq")
# BiocManager::install("edgeR")
# BiocManager::install("EDASeq")

library (SummarizedExperiment)
library (biomaRt)
library (vroom)
library (NOISeq)
library (edgeR)

# Leer data --------------------------------------------------------------------

healthy <- vroom (file = 'FPKM_matrix_healthy.csv')
AD <- vroom (file = "AD.csv")
ExpeGroups <- vroom (file = "groups.csv")


# Crear mi ExpressionData ------------------------------------------------------

expre <- as.matrix(cbind(healthy[,2:45], AD [,2:135]))
rownames(expre) <- AD$gene_id  #aun se le tiene que agregar las anotaciones

# Anotaciones biologicas adicionales / GC content, GeneBiotype, HGNCinfo, Length

mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values=rownames(expre),mart=mart)

myannot$length=abs(myannot$end_position-myannot$start_position)

# Filtrar los transcritos sin anotaciones

myannot=myannot[myannot$gene_biotype=="protein_coding"&
                  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
ExpreData =expre[rownames(expre)%in%myannot$ensembl_gene_id,] #ExpressionData
dim(ExpreData)
# 18848   178

# Convertir a NOISeq object ----------------------------------------------------

noiseqData = NOISeq::readData(data = ExpreData,
                      factor = as.data.frame(ExpeGroups),
                      length = myannot[,c(1,8)],
                      biotype = myannot[,c(1,3)],
                      gc = myannot[,1:2])

########  CHECK BIASES ---------------------------------------------------------

# 1. Counts distribution per biotype 
mycountsbio = dat(noiseqData, type = "countsbio", factor = "group")
explo.plot(mycountsbio, plottype = "barplot")
# 2. Check for low counts genes
hist(rowMeans(cpm(ExpreData,log=T)),ylab="genes",
     xlab="mean of log CPM",col="gray")
# 3. Check RNA composition bias
mycd = dat(noiseqData, type = "cd", norm = TRUE) #slooooow
# "Diagnostic test: FAILED. Normalization is required to correct this bias."
# "Confidence intervals for median of M:"
# 0.01%                   99.99%   
explo.plot(mycd,samples=sample(1:ncol(ExpreData),10))
# 4. Check for GC and length bias 
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias",
                   factor = "group")
explo.plot(myGCcontent, samples = NULL, toplot = "global")
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias",
                 factor = "group")
explo.plot(mylenBias, samples = NULL, toplot = "global")
#5) check for batch effect
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = F)
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "group")

######## SOLVE BIASES ----------------------------------------------------------

library(EDASeq)
# 1. Filter low count genes/ CPM Method
# CPM =(counts/fragments sequenced)*one million.

ExpreDataFiltered = filtered.data(ExpreData, factor = "group",
                                    norm = TRUE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#14959 features are to be kept for differential expression analysis with filtering method 1

myannot=myannot[myannot$ensembl_gene_id%in%rownames(ExpreDataFiltered),]

#correcting rownames order
myannot <- myannot[order(myannot$ensembl_gene_id), ]
dim(myannot)
# [1] 14959     8
ExpreDataFiltered <- ExpreDataFiltered[match(myannot$ensembl_gene_id, rownames(ExpreDataFiltered)), ] 
dim(ExpreDataFiltered)
# [1] 14959   178
identical(row.names(as.matrix(ExpreDataFiltered)),row.names(data.frame(myannot,row.names=myannot$ensembl_gene_id)))

#all names must match
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(ExpreDataFiltered),
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(ExpeGroups,row.names=ExpeGroups$sample))

#order for less bias
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content", which = "full") #corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full") #corrects length bias 

# UQUA normalization
fullfullUQUA <-NOISeq::uqua(normCounts(lFull), long = 1000, lc = 0, k = 0)
norm.counts <- betweenLaneNormalization(normCounts(lFull), which = "median", offset = FALSE) # ------ here------
noiseqData = NOISeq::readData(data = fullfullUQUA, factors=as.data.frame(ExpeGroups))

#cd has to preceed ARSyN or won't work
mycd=NOISeq::dat(noiseqData,type="cd",norm=TRUE)
# [1] "Diagnostic test: PASSED."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
##PASSED KEI 
# 177

####### SOLVE BATCH EFFECT -----------------------------------------------------
myPCA = dat(noiseqData, type = "PCA", norm = T, logtransf = F)
explo.plot(myPCA, samples = c(1,2), plottype = "scores",
           factor = "group")
ffTMMARSyn=ARSyNseq(noiseqData, factor = "group", batch = F,
                    norm = "n",  logtransf = T)
myPCA = dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = T)
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
           factor = "group")

######## FINAL QUALITY CHECK ---------------------------------------------------
expr_mat <- exprs (ffTMMARSyn) #matriz de expresion buena 
saveRDS(expr_mat, "ADUQUARSyn.rds")

noiseqData = NOISeq::readData(data = expr_mat,
                              factor = as.data.frame(ExpeGroups),
                              length = myannot[,c(1,8)],
                              biotype = myannot[,c(1,3)],
                              gc = myannot[,1:2])

mycountsbio = dat(noiseqData, type = "countsbio", factor = "group",
                  norm=T)
explo.plot(mycountsbio, plottype = "boxplot", samples=1:2)
explo.plot(mycountsbio, toplot = 2, samples = 1:2, plottype = "boxplot")
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias", 
                   factor = "group",norm=T)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(myGCcontent, samples = x))
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", 
                 factor = "group",norm=T)
par(mfrow=c(1,2))
sapply(1:2,function(x) explo.plot(mylenBias, samples = x))

# END --------------------------------------------------------------------------
