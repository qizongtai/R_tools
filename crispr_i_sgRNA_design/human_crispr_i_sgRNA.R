#######################################################
#design crispr_i human sgRNAs from Dolcetto lab and order oligo from IDT
#https://www.addgene.org/pooled-library/broadgpp-human-crispri-dolcetto/
#Choose "Human Dolcetto Set A Target Genes"
#https://media.addgene.org/cms/filer_public/1c/59/1c59fe51-ef6e-44bf-963a-598edadcd66f/broadgpp-dolcetto-targets-seta.txt

#step1: get the crispr_i guide sequences
library(readxl)
path <- "C:/Users/qizon/Dropbox/12-R working dir/6_snai2"
setwd(path)
allguide <- read.table(file = "1_data/broadgpp-dolcetto-targets-seta.txt", sep="\t", header = T)
gene <- read_xlsx(path = "1_data/scc47_and_jhu6_dw_set1.xlsx", col_names = T, sheet = 1)
gene.guide <- merge(x=allguide, y=gene[1], by.x=2, by.y=1)
#add two other columns
gene.guide.test <- data.frame(gene.guide, "Scale" = rep("25nm", nrow(gene.guide)), "Purification" = rep("STD", nrow(gene.guide)))
#gene.guide[,c("Scale", "Purification")] <- data.frame(rep("25nm", nrow(gene.guide)), rep("STD", nrow(gene.guide)))
#order two guides so delete the third one 
gene.guide.2 <- gene.guide[-seq(from = 3, to = 60, by = 3),]
write.table(gene.guide.2, file = "3_output/dw_gene_guides.xls", sep = "\t", quote = F, col.names = T, row.names = F)
#write.xlsx(gene.guide.2, sheetName = "3_output/dw_gene_guides.xlsx", col.names = T, row.names = F)

#step2: make annealing primers, transfer to the IDT excel file
library(Biostrings)
anneal_fwd <- function(x) {return(paste0("CACCG",x))}
anneal_rev <- function(x) {return(paste0("AAAC",reverseComplement(DNAString(x)),"C"))}
l <- list()
for (i in 1:nrow(gene.guide.2)) {
  n <- ifelse(i %% 2 == 1, 1, 2)
  l[paste0(gene.guide.2[i,1], "_", n, "_fwd")] <- anneal_fwd(gene.guide.2[i,2])
  l[paste0(gene.guide.2[i,1], "_", n, "_rev")] <- anneal_rev(gene.guide.2[i,2])
}
m <- as.matrix(l)
write.table(m, file = "3_output/dw_gene_guides_primers.xls", col.names = F, row.names = T, quote = F, sep = "\t")

#######################################################
#design crispr_i human sgRNAs from Dolcetto lab and order oligo from IDT
#https://www.addgene.org/pooled-library/broadgpp-human-crispri-dolcetto/
#Choose "Human Dolcetto Set A Target Genes"
#https://media.addgene.org/cms/filer_public/1c/59/1c59fe51-ef6e-44bf-963a-598edadcd66f/broadgpp-dolcetto-targets-seta.txt

#step1: get the crispr guide sequences
library(readxl)
library(data.table)
path <- "C:/Users/qizon/Dropbox/12-R working dir/6_snai2"
setwd(path)
allguide <- read.table(file = "1_data/broadgpp-dolcetto-targets-seta.txt", sep="\t", header = T)
emt.tf <- data.frame(gene=c("SNAI1", "TWIST1", "TWIST2", "ZEB1", "ZEB2"))
gene.guide <- merge(x=emt.tf, y=allguide, by.x = "gene", by.y=2)

#step2: make annealing primers, transfer to the IDT excel file
library(Biostrings)
anneal_fwd <- function(x) {return(paste0("CACCG",x))}
anneal_rev <- function(x) {return(paste0("AAAC",reverseComplement(DNAString(x)),"C"))}

#option1: get just two out of three
gene.guide.2 <- gene.guide[-seq(from = 3, to = 60, by = 3),]
l <- list()
for (i in 1:nrow(gene.guide.2)) {
  n <- ifelse(i %% 2 == 1, 1, 2)
  l[paste0(gene.guide.2[i,1], "_", n, "_fwd")] <- anneal_fwd(gene.guide.2[i,2])
  l[paste0(gene.guide.2[i,1], "_", n, "_rev")] <- anneal_rev(gene.guide.2[i,2])
}
m <- as.matrix(l)
write.table(m, file = "3_output/emt_gene_2guides_primers_.xls", col.names = F, row.names = T, quote = F, sep = "\t")

#option2: get all three primers
l <- list()
for (i in 1:nrow(gene.guide)) {
  n <- list("1"=1, "2"=2, "0"=3)
  l[paste0(gene.guide[i,1], "_", n[[as.character(i %% 3)]], "_fwd")] <- anneal_fwd(gene.guide[i,2])
  l[paste0(gene.guide[i,1], "_", n[[as.character(i %% 3)]], "_rev")] <- anneal_rev(gene.guide[i,2])
}
df <- as.data.frame(l)
df.t <- transpose(df)
rownames(df.t) <- colnames(df)
write.table(df.t, file = "3_output/emt_gene_3guides_primers.xls", col.names = F, row.names = T, quote = F, sep = "\t")
