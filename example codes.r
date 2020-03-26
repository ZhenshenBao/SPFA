library(KEGGgraph)
library(igraph)
library(limma)
library(ggplot)

set = "GSE4183"
data(list = set, package = "KEGGdzPathwaysGEO")
sets <- get(set)
#write a function to extract required info into a list
getdataaslist = function(x) {
  x = get(x)
  exp = experimentData(x)
  dataset = exp@name
  disease = notes(exp)$disease
  dat.m = exprs(x)
  ano = pData(x)
  design = notes(exp)$design
  annotation = paste(x@annotation, ".db", sep = "")
  targetGeneSets = notes(exp)$targetGeneSets
  list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
  names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
                  "targetGeneSets")
  return(list)
}
dlist = getdataaslist(set)

block = dlist$ano$Block
Block = block
block = factor(Block)
group = dlist$ano$Group
G = factor(group)
force(G)
force(block)
paired = dlist$design == "Paired"
esetm = dlist$dat.m
annotation = dlist$annotation

if (paired) {
  stopifnot(length(block) == length(group))
  stopifnot(all(table(block) == 2))
}

aT1 = filteranot(esetm, group, paired, block, annotation)
esetm = esetm[rownames(esetm) %in% aT1$ID, ]
rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]

topSigNum = dim(esetm)[1]

if (paired) {
  design <- model.matrix(~0 + G + block)
  colnames(design) <- substr(colnames(design), 2, 100)
}
if (!paired) {
  design <- model.matrix(~0 + G)
  colnames(design) <- levels(G)
}
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts = "d-c", levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1 <- topTable(fit2, coef = 1, number = topSigNum)
aT1$ID <- rownames(aT1)
head(aT1)

all <-  aT1$ID
all[1:10]

de <- aT1[which(aT1$adj.P.Val<0.05),]
de <- de$ID
if(length(de) < 200){
  top <- aT1[which(aT1$P.Value<0.05),]
  top <- top[which(top$logFC>1.5),]
  de <- top$ID
}
if(length(de) < 200){
  top <- aT1[1:(nrow(aT1)*0.01),]
  de <- top$ID
}
de[1:10]

normalsamples <- dlist$ano[which(dlist$ano$Group=="c"),1]
normal <- esetm[,normalsamples]
cancersamples <- dlist$ano[which(dlist$ano$Group=="d"),1]
cancer <- esetm[,cancersamples]

normal <- t(normal)
colnames(normal) <- as.character(colnames(normal))
normal <- as.data.frame(normal)
head(normal)

cancer <- t(cancer)
colnames(cancer) <- as.character(colnames(cancer))
cancer <- as.data.frame(cancer)
head(cancer)

res <- SPFA(paths=paths,mdir=mdir,de=de,all=all,normal=normal,cancer=cancer)
head(res)

res_per <- SPFA_per(mdir=mdir,pathfilename="hsa04151.xml",normal=normal,cancer=cancer)
res_per
