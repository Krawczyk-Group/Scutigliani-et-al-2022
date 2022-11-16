if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("GEOquery", force = TRUE)
BiocManager::install("biomaRt")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomeInfoDb")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("hgu133plus2.db")
library(GEOquery)
library(Biobase)
library("biomaRt")
library("dplyr")
library("AnnotationDbi")
library("org.Hs.eg.db")
#library("illuminaHumanv4.db")
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)



#####    Tabuchi et al., 2008   #############

#Read the dataset
gset <- getGEO("GSE10043", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)

#Define the vectors where the data will be saved
N<-dim(ex)[1]
name1<-rep(N)
r1<-rep(0, N)

#Loop to save the genes with ratios

for (i in 1:N){
  name1[i]<-rownames(exprs(gset))[i]
  r1[i]<-(((ex[i,3]/ex[i,1])+(ex[i,4]/ex[i,2]))/2)

}

#Change ensamble

ch <- getBM(mart=mart, attributes=c("affy_hg_u133_plus_2","ensembl_gene_id", "gene_biotype", "external_gene_name"), filter="affy_hg_u133_plus_2", values <- name1 , uniqueRows=TRUE)

#Vectors with the old notation (not) and the gen symbols (symb)
not<-ch[,1]
eid<-ch[,2]
symbo<-ch[,4]

#Order the gene symbol vector so it is in the same order as the vector with the gene regulations
N<-length(eid)[1]
reg1<-rep(0, length(eid))
eid1<-rep(0, length(eid))
symbo1<-rep(0, length(eid))
affy1<-rep(0, length(eid))

N<-length(name1)
j<-1
for(i in 1:N){
  n<-length(which(not==name1[i]))
  index<-which(not==name1[i])
  if(n>0){
    for(k in 1:n){
      eid1[j]<-eid[index[k]]
      symbo1[j]<-symbo[index[k]]
      reg1[j]<-r1[i]
      affy1[j]<-not[index[k]]
      j<-j+1
    }
  }else{
    if(n!=0){
      eid1[j]<-eid[index]
      symbo1[j]<-symbo[index]
      affy1[j]<-not[index]
      reg1[j]<-r1[i]
      j<-j+1
    }
  }
}

## Removing Non coding DNA
reg1<-reg1[!duplicated(affy1)]
eid1<-eid1[!duplicated(affy1)]
symbo1<-symbo1[!duplicated(affy1)]
affy1<-affy1[!duplicated(affy1)]
affy1<-affy1[symbo1!=""]
reg1<-reg1[symbo1!=""]
eid1<-eid1[symbo1!=""]
symbo1<-symbo1[symbo1!=""]

#Create a single vector with all the gene symbols which are present in the datasets (without repetitions)
name<-c(eid1)
name<-name[!is.na(name)]

#Define the vectors in which the final results will be saved
N<-length(name)

#Create a matrix to count the number of papers(column 1)and number of datasets(column 2) in which a gen is present as well how many times it is upregulated(column3), downregulated(column4) and noregulated(column5)

x <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x) <- c("Gene_Symbol","Tabuchi2008","log")

#rownames(repetition) <- name
for(i in 1:N){
  x[i,1]<-symbo1[i]
  x[i,2]<-reg1[i]
  x[i,3]<-(log10(reg1[i]))^2
  
}

x1.frame <- as.data.frame(x)
x1.frame <- x1.frame[with(x1.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x1.frame = subset(x1.frame, select = -c(log) )

####   Tabuchi et al., 2011   ################

#Read the dataset
gset <- getGEO("GSE23405", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)


#Define the vectors where the data will be saved
N<-dim(ex)[1]
name2<-rep(0, N)
r21<-rep(0, N)
r22<-rep(0, N)
r23<-rep(0, N)
r24<-rep(0, N)
r25<-rep(0, N)
r26<-rep(0, N)
r27<-rep(0, N)
r28<-rep(0, N)


#Loop to save the genes with ratios

for (i in 1:N){
  name2[i]<-rownames(exprs(gset))[i]
  r21[i]<-(ex[i,2]/ex[i,1])
  r22[i]<-(ex[i,3]/ex[i,1])
  r23[i]<-(ex[i,4]/ex[i,1])
  r24[i]<-(ex[i,5]/ex[i,1])
  r25[i]<-(ex[i,6]/ex[i,1])
  r26[i]<-(ex[i,7]/ex[i,1])
  r27[i]<-(ex[i,8]/ex[i,1])
  r28[i]<-(ex[i,9]/ex[i,1])
}

#Change ensamble

ch <- getBM(mart=mart, attributes=c("affy_hg_u133_plus_2","ensembl_gene_id", "gene_biotype", "external_gene_name"), filter="affy_hg_u133_plus_2", values <- name2 , uniqueRows=TRUE)



#Vectors with the old notation (not) and the gen symbols (symb)
not<-ch[,1]
eid<-ch[,2]
symbo<-ch[,4]

#Order the gene symbol vector so it is in the same order as the vector with the gene regulations
N<-length(eid)[1]
reg21<-rep(0, length(eid))
reg22<-rep(0, length(eid))
reg23<-rep(0, length(eid))
reg24<-rep(0, length(eid))
reg25<-rep(0, length(eid))
reg26<-rep(0, length(eid))
reg27<-rep(0, length(eid))
reg28<-rep(0, length(eid))
eid2<-rep(0, length(eid))
symbo2<-rep(0, length(eid))
affy2<-rep(0, length(eid))



N<-length(name2)
j<-1
for(i in 1:N){
  n<-length(which(not==name2[i]))
  index<-which(not==name2[i])
  if(n>0){
    for(k in 1:n){
      eid2[j]<-eid[index[k]]
      symbo2[j]<-symbo[index[k]]
      affy2[j]<-not[index[k]]
      reg21[j]<-r21[i]
      reg22[j]<-r22[i]
      reg23[j]<-r23[i]
      reg24[j]<-r24[i]
      reg25[j]<-r25[i]
      reg26[j]<-r26[i]
      reg27[j]<-r27[i]
      reg28[j]<-r28[i]
      j<-j+1
    }
  }else{
    if(n!=0){
      eid2[j]<-eid[index]
      symbo2[j]<-symbo[index]
      affy2[j]<-not[index]
      reg21[j]<-r21[i]
      reg22[j]<-r22[i]
      reg23[j]<-r23[i]
      reg24[j]<-r24[i]
      reg25[j]<-r25[i]
      reg26[j]<-r26[i]
      reg27[j]<-r27[i]
      reg28[j]<-r28[i]
      j<-j+1
    }
  }
}

## Removing Non coding DNA
reg21<-reg21[!duplicated(affy2)]
reg22<-reg22[!duplicated(affy2)]
reg23<-reg23[!duplicated(affy2)]
reg24<-reg24[!duplicated(affy2)]
reg25<-reg25[!duplicated(affy2)]
reg26<-reg26[!duplicated(affy2)]
reg27<-reg27[!duplicated(affy2)]
reg28<-reg28[!duplicated(affy2)]
eid2<-eid2[!duplicated(affy2)]
symbo2<-symbo2[!duplicated(affy2)]
affy2<-affy2[!duplicated(affy2)]
affy2<-affy2[symbo2!=""]
reg21<-reg21[symbo2!=""]
reg22<-reg22[symbo2!=""]
reg23<-reg23[symbo2!=""]
reg24<-reg24[symbo2!=""]
reg25<-reg25[symbo2!=""]
reg26<-reg26[symbo2!=""]
reg27<-reg27[symbo2!=""]
reg28<-reg28[symbo2!=""]
eid2<-eid2[symbo2!=""]
symbo2<-symbo2[symbo2!=""]

#Create a single vector with all the gene symbols which are present in the datasets (without repetitions)
name<-c(eid2)
name<-name[!is.na(name)]

#Define the vectors in which the final results will be saved
N<-length(name)

#Create a matrix to count the number of papers(column 1)and number of datasets(column 2) in which a gen is present as well how many times it is upregulated(column3), downregulated(column4) and noregulated(column5)

x21 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x21) <- c("Gene_Symbol","Tabuchi2011a","log")

for(i in 1:N){
  x21[i,1]<-symbo2[i]
  x21[i,2]<-reg21[i]
  x21[i,3]<-(log10(reg21[i]))^2
  
}

x21.frame <- as.data.frame(x21)
x21.frame <- x21.frame[with(x21.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x21.frame = subset(x21.frame, select = -c(log) )

x22 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x22) <- c("Gene_Symbol","Tabuchi2011b","log")

for(i in 1:N){
  x22[i,1]<-symbo2[i]
  x22[i,2]<-reg22[i]
  x22[i,3]<-(log10(reg22[i]))^2
  
}

x22.frame <- as.data.frame(x22)
x22.frame <- x22.frame[with(x22.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x22.frame = subset(x22.frame, select = -c(log) )

x23 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x23) <- c("Gene_Symbol","Tabuchi2011c","log")

for(i in 1:N){
  x23[i,1]<-symbo2[i]
  x23[i,2]<-reg23[i]
  x23[i,3]<-(log10(reg23[i]))^2
  
}

x23.frame <- as.data.frame(x23)
x23.frame <- x23.frame[with(x23.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x23.frame = subset(x23.frame, select = -c(log) )


x24 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x24) <- c("Gene_Symbol","Tabuchi2011d","log")

for(i in 1:N){
  x24[i,1]<-symbo2[i]
  x24[i,2]<-reg24[i]
  x24[i,3]<-(log10(reg24[i]))^2
  
}

x24.frame <- as.data.frame(x24)
x24.frame <- x24.frame[with(x24.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x24.frame = subset(x24.frame, select = -c(log) )

x25 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x25) <- c("Gene_Symbol","Tabuchi2011e","log")

for(i in 1:N){
  x25[i,1]<-symbo2[i]
  x25[i,2]<-reg25[i]
  x25[i,3]<-(log10(reg25[i]))^2
  
}

x25.frame <- as.data.frame(x25)
x25.frame <- x25.frame[with(x25.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x25.frame = subset(x25.frame, select = -c(log) )

x26 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x26) <- c("Gene_Symbol","Tabuchi2011f","log")

for(i in 1:N){
  x26[i,1]<-symbo2[i]
  x26[i,2]<-reg26[i]
  x26[i,3]<-(log10(reg26[i]))^2
  
}

x26.frame <- as.data.frame(x26)
x26.frame <- x26.frame[with(x26.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x26.frame = subset(x26.frame, select = -c(log) )

x27 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x27) <- c("Gene_Symbol","Tabuchi2011g","log")

for(i in 1:N){
  x27[i,1]<-symbo2[i]
  x27[i,2]<-reg27[i]
  x27[i,3]<-(log10(reg27[i]))^2
  
}

x27.frame <- as.data.frame(x27)
x27.frame <- x27.frame[with(x27.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x27.frame = subset(x27.frame, select = -c(log) )
x28 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x28) <- c("Gene_Symbol","Tabuchi2011h","log")

for(i in 1:N){
  x28[i,1]<-symbo2[i]
  x28[i,2]<-reg28[i]
  x28[i,3]<-(log10(reg28[i]))^2
  
}

x28.frame <- as.data.frame(x28)
x28.frame <- x28.frame[with(x28.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x28.frame = subset(x28.frame, select = -c(log) )


####  Court, et al 2017 #####


#Read the dataset
gset <- getGEO("GSE92990", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)

#Define the vectors where the data will be saved
N<-dim(ex)[1]
name3<-rep(0, N)
r3<-rep(0, N)

#Loop to save the genes with ratios

for (i in 1:N){
  name3[i]<-rownames(exprs(gset))[i]
  r3[i]<-(((ex[i,1]/ex[i,4])+(ex[i,2]/ex[i,5])+(ex[i,3]/ex[i,6]))/3)

}

#Change the name of the genes to their gene symbol
ch<-data.frame(Gene=unlist(mget(x = name3,envir = illuminaHumanv4SYMBOL)))
#ch<???na.omit(ch) 

#Delete the cases in which the gene symbol could not be found 
not<-row.names(ch)
symb3<- as.vector(unlist(ch$Gene))

#Order the gene symbol vector so it is in the same order as the vector with the gene regulations
N<-length(symb3)[1]
reg3<-rep(0, length(symb3))
symbo3<-rep(0, length(symb3))


N<-length(name3)
j<-1
for(i in 1:N){
  n<-length(which(not==name3[i]))
  index<-which(not==name3[i])
  if(n>0){
    for(k in 1:n){
      symbo3[j]<-symb3[index[k]]
      reg3[j]<-r3[i]
      j<-j+1
    }
  }else{
    if(n!=0){
      symbo3[j]<-symb3[index]
      reg3[j]<-r3[i]
      j<-j+1
    }
  }
}

## Removing Non coding DNA
symbo3<-symbo3[symbo3!=""]
reg3<-reg3[symbo3!=""]

#Create a single vector with all the gene symbols which are present in the datasets (without repetitions)
name<-c(symbo3)
name<-name[!is.na(name)]

#Define the vectors in which the final results will be saved
N<-length(name)

#Create a matrix to count the number of papers(column 1)and number of datasets(column 2) in which a gen is present as well how many times it is upregulated(column3), downregulated(column4) and noregulated(column5)

x <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x) <- c("Gene_Symbol","Court2017","log")

#rownames(repetition) <- name
for(i in 1:N){
  x[i,1]<-symbo3[i]
  x[i,2]<-reg3[i]
  x[i,3]<-(log10(reg3[i]))^2
  
}

x3.frame <- as.data.frame(x)
x3.frame <- x3.frame[with(x3.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x3.frame = subset(x3.frame, select = -c(log) )
x3.frame <- na.omit(x3.frame )




####   Amaya et al., 2014   #####

#Read the dataset
gset <- getGEO("GSE48398", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)


#Define the vectors where the data will be saved
N<-dim(ex)[1]
name4<-rep(0, N)
r41<-rep(0, N)
r42<-rep(0, N)
r43<-rep(0, N)


#Loop to save the genes with ratios

for (i in 1:N){
  name4[i]<-rownames(exprs(gset))[i]
  r41[i]<-(mean(ex[i,14:16])/mean(ex[i,7:13]))
  r42[i]<-(mean(ex[i,23:25])/mean(ex[i,17:22]))
  r43[i]<-(mean(ex[i,31:36])/mean(ex[i,25:30]))
}

#Change the name of the genes to their gene symbol
ch<-data.frame(Gene=unlist(mget(x = name4,envir = illuminaHumanv4SYMBOL)))

not<-row.names(ch)
symbo<- as.vector(unlist(ch$Gene))



#Order the gene symbol vector so it is in the same order as the vector with the gene regulations
N<-length(symbo)[1]
reg41<-rep(0, length(symbo))
reg42<-rep(0, length(symbo))
reg43<-rep(0, length(symbo))
symbo4<-rep(0, length(symbo))



N<-length(name4)
j<-1
for(i in 1:N){
  n<-length(which(not==name4[i]))
  index<-which(not==name4[i])
  if(n>0){
    for(k in 1:n){
      symbo4[j]<-symbo[index[k]]
      reg41[j]<-r41[i]
      reg42[j]<-r42[i]
      reg43[j]<-r43[i]
      j<-j+1
    }
  }else{
    if(n!=0){
      symbo4[j]<-symbo[index]
      reg41[j]<-r41[i]
      reg42[j]<-r42[i]
      reg43[j]<-r43[i]
      j<-j+1
    }
  }
}

## Removing Non coding DNA
symbo4<-symbo4[symbo4!=""]
reg41<-reg41[symbo4!=""]
reg42<-reg42[symbo4!=""]
reg43<-reg43[symbo4!=""]

#Create a single vector with all the gene symbols which are present in the datasets (without repetitions)
name<-c(symbo4)
#name<-name[!is.na(name)]

#Define the vectors in which the final results will be saved
N<-length(name)

#Create a matrix to count the number of papers(column 1)and number of datasets(column 2) in which a gen is present as well how many times it is upregulated(column3), downregulated(column4) and noregulated(column5)

x41 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x41) <- c("Gene_Symbol","Amaya2014a","log")

#rownames(repetition) <- name
for(i in 1:N){
  x41[i,1]<-symbo4[i]
  x41[i,2]<-reg41[i]
  x41[i,3]<-(log10(reg41[i]))^2
  
}

x41.frame <- as.data.frame(x41)
x41.frame <- x41.frame[with(x41.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x41.frame = subset(x41.frame, select = -c(log) )
x41.frame <- na.omit(x41.frame)

x42 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x42) <- c("Gene_Symbol","Amaya2014a","log")

#rownames(repetition) <- name
for(i in 1:N){
  x42[i,1]<-symbo4[i]
  x42[i,2]<-reg42[i]
  x42[i,3]<-(log10(reg42[i]))^2
  
}

x42.frame <- as.data.frame(x42)
x42.frame <- x42.frame[with(x42.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x42.frame = subset(x42.frame, select = -c(log) )
x42.frame <- na.omit(x42.frame)

x43 <- matrix(rep(0,3*N),nrow = N, ncol = 3)
colnames(x43) <- c("Gene_Symbol","Amaya2014a","log")

#rownames(repetition) <- name
for(i in 1:N){
  x43[i,1]<-symbo4[i]
  x43[i,2]<-reg43[i]
  x43[i,3]<-(log10(reg43[i]))^2
  
}

x43.frame <- as.data.frame(x43)
x43.frame <- x43.frame[with(x43.frame, ave(log, Gene_Symbol, FUN=max)==log),]
x43.frame = subset(x43.frame, select = -c(log) )
x43.frame <- na.omit(x43.frame)

library(xlsx) 
write.csv(x1.frame,"Tabuchi2008.csv", row.names = FALSE)
write.csv(x21.frame,"Tabuchi2011a.csv", row.names = FALSE)
write.csv(x22.frame,"Tabuchi2011b.csv", row.names = FALSE)
write.csv(x23.frame,"Tabuchi2011c.csv", row.names = FALSE)
write.csv(x24.frame,"Tabuchi2011d.csv", row.names = FALSE)
write.csv(x25.frame,"Tabuchi2011e.csv", row.names = FALSE)
write.csv(x26.frame,"Tabuchi2011f.csv", row.names = FALSE)
write.csv(x27.frame,"Tabuchi2011g.csv", row.names = FALSE)
write.csv(x28.frame,"Tabuchi2011h.csv", row.names = FALSE)
write.csv(x3.frame,"Court2017.csv", row.names = FALSE)
write.csv(x41.frame,"Amaya2014a.csv", row.names = FALSE)
write.csv(x42.frame,"Amaya2014b.csv", row.names = FALSE)
write.csv(x43.frame,"Amaya2014c.csv", row.names = FALSE)

