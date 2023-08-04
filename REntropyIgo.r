## load libraries
library(SLICE)   #loading the SLICE library
library(anndata)
library(Biobase)

## define inputs
data.name = "paul_15"
species = "mouse"
celltype = "paul15_clusters"

## Build ExpressionSet

# phenoData:
tmp <- read.csv("adataobs.csv", row.names = 1) #import obs
tmp$Cell = rownames(tmp)
pdata <- AnnotatedDataFrame(tmp)

print(tmp)
# featureData:
tmp <- read.csv("adatavar.csv", row.names = 1) #import var
tmp$SYMBOL = rownames(tmp)
fdata <- AnnotatedDataFrame(tmp)

print(tmp)
# expression data:
tmp <- read.table("gene_expression.csv", row.names = 1, header = TRUE, sep = ",") #import X
m <- as.matrix(tmp)
m = t(m) 
row.names(m) = fdata@data$SYMBOL

## create ExpressionSet object:
eset <- new("ExpressionSet", exprs = m, phenoData = pdata, featureData = fdata)


## Run SLICE
es = eset

context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="")

ercc.genes <- grep("^ERCC-", rownames(fData(es)), value = TRUE)
rb.genes <- grep("^Rpl|^Rps|^Mrpl|^Mrps", rownames(fData(es)), value = TRUE)
es <- es[which(!(rownames(fData(es)) %in% c(rb.genes, ercc.genes))), ]

sc <- construct(exprmatrix=as.data.frame(exprs(es)), 
                cellidentity=factor(pData(es)[celltype]),
                projname=context_str)

data(mm_kappasim)

sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=100,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=201602)                    # set the random seed to reproduce the results in the paper

I_go = sc@entropies
colnames(I_go) = "I_go"

## Export results
write.csv(I_go, "data/Rigo.csv",row.names=FALSE)