#################################################################################
## Proyecto: Mecanismos de desregulación del metabolismo energético asociados a 
##           cáncer de mama

################################################################################
## EXPRESIÓN DIFERENCIAL
## Autor: Cristóbal Fresno
## Fecha: 2019/02/05
## Modificado por : Erandi Serrano
## https://github.com/CSB-IG/tcgarnaseqbc/blob/master/DifGenes.R
#################################################################################
#options(width=120)

# Instalar las bibliotecas necesarias

# NOIseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#    BiocManager::install("NOISeq", version = "3.8")

# EDAseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#      BiocManager::install("EDASeq", version = "3.8")
#      BiocManager::install("Glimma", version = "3.8")

library("NOISeq")
library("EDASeq")
library("Glimma")
      
CompleteMatrix <-  read.delim("CompleteMatrix.tsv", stringsAsFactors = FALSE)
CompleteMatrix <-  CompleteMatrix[-1,]
row.names(CompleteMatrix) <- CompleteMatrix[,1]
CompleteMatrix <- CompleteMatrix[,-1]
CompleteMatrix <- log2(CompleteMatrix)
#Insatalar limma

#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("limma", version = "3.8")
##Haciendo el análisis con limma
library("limma")

table(regexpr(pattern="N", names(CompleteMatrix))>0)
# FALSE  TRUE 
#   735   113 
muestras <- data.frame(Sanos = regexpr(pattern="N", names(CompleteMatrix))>0)
muestras$Grupo <- "S"
muestras$Grupo[!muestras$Sanos]  <- "T"
design<-model.matrix(~1+Grupo, data=muestras)
head(design)
#     (Intercept) GrupoTCompleteMatrix.tsv
# 1             1      0
# 2             1      0
# 847           1      1
# 848           1      1
# y= mu + Tumor + Error

fit <- lmFit(CompleteMatrix, design)
head(fit$coefficients)
#           (Intercept)     GrupoT
# A4GALT     10.167763 -0.6019043
# AACS       11.536646 -0.1389764
# AADAT       7.335554 -0.3976854
# AARS       12.579470  0.3040321
# AARS2      10.471358  0.1187261
# AASDHPPT   11.399976  0.1607934

##Enfermos-Sanos están en GroupT
fit2 <- eBayes(fit)
fit2$fdr<-apply(fit2$"p.value", 2, p.adjust, method="fdr")

library("ggplot2")
library("reshape2")

p<-ggplot(as.data.frame(fit$coefficients), aes(x=GrupoT))+geom_density()
pdf(file="LFCDensity.pdf")
pdegCount
dev.off()

nrow(CompleteMatrix)
# [1] 1298

##Buscando los genes diferenciales con ebayes---------------------------------
alphas<-c(0.05, 10^(-2:-10), 10^(seq(-20, -100, by=-10)))
degCount<-sapply(alphas, function(alpha){
  table(fit2$fdr[,"GrupoT"]<alpha)
})
colnames(degCount)<-alphas
degCount<-as.data.frame(t(degCount))
degCount$alpha<-alphas
degCount$Noise<-alphas*nrow(CompleteMatrix)
degCount

write.table(fit2, file = "20190207_DifExpMatrix.tsv", sep = "\t", row.names = F)
write.table(CompleteMatrix, file = "20190207_CompleteMatrix_log2.tsv", sep = "\t", row.names = F)

#Volver a pegar columna de los nombres de genes en matriz de log2
#TREAT
lfc<-seq(0, 3, by=0.5)
B <- 5
fitTreat<-lapply(lfc, function(x){
  ##Ajustando el modelo
  aux<-treat(fit, lfc=x)
  aux$fdr<-apply(aux$"p.value", 2, p.adjust, method="fdr")
  ##Estadistico B
  p<-1-aux$fdr[, "GrupoT"] ##Probabilidad de salir diferencialºº
  aux$B<-log(p/(1-p))
  ##Diferenciales
  aux$dif<-data.frame(
    alpha=alphas,
    lfc=rep(x, length(alphas)),
    Dif=t(sapply(alphas, function(alpha){
      table(factor(aux$fdr[,"GrupoT"]<alpha & (aux$B>5), levels=c(TRUE, FALSE)))
    })),
    Noise=alphas*nrow(CompleteMatrix)
  )
  
  return(aux)
})
names(fitTreat)<-paste("lfc", lfc, sep="")
difTreat<-do.call(rbind, lapply(fitTreat, function(x){x$dif}))
head(difTreat)

#Glimma plot
fitTreat$lfc1.5$deg<-fitTreat$lfc1.5$fdr[,"GrupoT", drop=FALSE]<0.01 & (fitTreat$lfc1.5$B>5)
table(fitTreat$lfc1.5$deg)
dt <- decideTests(fitTreat$lfc1.5)
#La ocupamos despues para pintar las redes y por eso la guardamos...
write.table(dt, file = "dt.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = FALSE )
summary(dt)
glMDPlot(fitTreat$lfc1.5, counts = CompleteMatrix, groups = muestras$Grupo, status = dt)

#Salida
genesFULL<-cbind(
  Coef=fitTreat$lfc1.5$coefficients,
  Diff=fitTreat$lfc1.5$deg[, "GrupoT"], 
  "p.value"=do.call(cbind, lapply(fitTreat, function(x){x$"p.value"[, "GrupoT"]})),
  FDR=do.call(cbind, lapply(fitTreat, function(x){x$fdr[, "GrupoT"]})), 
  B=do.call(cbind, lapply(fitTreat, function(x){
    ##Estadistico B
    p<-1-x$fdr[, "GrupoT"] ##Probabilidad de salir diferencial
    B<-log(p/(1-p))
    return(B)
  })),
  Exp=fitTreat$lfc0.5$coefficients%*%t(unique(design)))

head(genesFULL)
genesFULL = genesFULL[,c("GrupoT", "lfc1.5")]
genesFULL = genesFULL[,c(fitTreat$lfc1.5, "lfc1.5")]
write.table(genesFULL, file = "genes_metabolic.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )
