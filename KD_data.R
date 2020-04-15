library(tools)

setwd("/Users/Denise/Documents/Master/IMBEI/data/14_APA_genes_whole_plus_dropouts")

temp = list.files(pattern=".txt")

# create empty matrices with a column for each condition
KD_matrix <- matrix(,nrow = 0, ncol = length(temp))
pValP_matrix <- matrix(,nrow = 0, ncol = length(temp))
pValD_matrix <- matrix(,nrow = 0, ncol = length(temp))

# create list with all conditions/KDs
all_conditions <- list()

for (file in temp) {
  filename = file_path_sans_ext(file)    
  str = strsplit(filename, '_')
  condition = str[[1]][2]
  all_conditions[[length(all_conditions)+1]] <- condition
}

# set colnames in matrices to conditions
colnames(KD_matrix) <- all_conditions
colnames(pValP_matrix) <- all_conditions
colnames(pValD_matrix) <- all_conditions

# read all condition files & create row in matrices for every gene, then fill in dir_idx and p-values for all affected genes
for (file in temp) {
  
  data <- read.table(file, header = TRUE, stringsAsFactors = TRUE)
  
  for (gene in data$Gene_symbol) {
    
    x <- character(length(temp))
    x <- NA
    
    if(!(gene %in% rownames(KD_matrix))) {
      KD_matrix <- rbind(KD_matrix,x)
      rownames(KD_matrix)[rownames(KD_matrix)=="x"] <- gene
      pValP_matrix <- rbind(pValP_matrix,x)
      rownames(pValP_matrix)[rownames(pValP_matrix)=="x"] <- gene
      pValD_matrix <- rbind(pValD_matrix,x)
      rownames(pValD_matrix)[rownames(pValD_matrix)=="x"] <- gene
    }
    
    filename = file_path_sans_ext(file)    
    str = strsplit(filename, '_')
    condition = str[[1]][2]
    
    #condition <- filename
    
    # extract dir index
    dir <- data[data$Gene_symbol == gene, "Dir_idx"]
    
    #exctract p-values
    pProx <- data[data$Gene_symbol == gene, "pvalue_prox"]
    pDist <- data[data$Gene_symbol == gene, "pvalue_dist"] 
    
    KD_matrix[gene, condition] <- dir
    pValP_matrix[gene, condition] <- pProx
    pValD_matrix[gene, condition] <- pDist
  }
  
}

# sort matrices by gene names (= rownames)
KD_matrix <- KD_matrix[order(rownames(KD_matrix)), ] 
pValP_matrix <- pValP_matrix[order(rownames(pValP_matrix)), ] 
pValD_matrix <- pValD_matrix[order(rownames(pValD_matrix)), ] 

# delete rows with all NAs (i.e. genes that are not affected in any condition)
KD_matrix_relevant <- KD_matrix[rowSums(is.na(KD_matrix)) != ncol(KD_matrix)-1,]
pValP_matrix_relevant <- pValP_matrix[rowSums(is.na(pValP_matrix)) != ncol(pValP_matrix)-1,]
pValD_matrix_relevant <- pValD_matrix[rowSums(is.na(pValD_matrix)) != ncol(pValD_matrix)-1,]

# turn zeros into NAs 
KD_matrix_relevant[KD_matrix_relevant == 0]  <- NA

setwd("/Users/Denise/Documents/Master/IMBEI/data/")

# save matrices in file
write.csv(KD_matrix_relevant, file = "KDmat.csv")
write.csv(pValP_matrix_relevant, file = "pProx.csv")
write.csv(pValD_matrix_relevant, file = "pDist.csv")


img.n <- readPNG("/Users/Denise/Documents/Master/IMBEI/data/plotCache/AARS_PAPD4_utr.png", native = TRUE)

plot(1:2, type= n )
rasterImage(img, 1.2, 1.27, 1.8, 1.73, interpolate=FALSE)

# list with genes
genenames <- rownames(KD_matrix_relevant)

# setwd("/Users/Denise/Documents/Master/IMBEI/data/")
# 
# # read KD data
# KDdata <- read.csv("KDmat.csv", dec = ".", sep = ",", header = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
# 
# KDmat <- as.matrix(KDdata)
# 
# str(KDmat)
# summary(KDmat)
# 
# # save conditions as colnames and genes as rownames
# colnames(KDmat) <- KDmat[1,]
# rownames(KDmat) <- KDmat[,1]
# colnames(KDmat)
# KDmat <- KDmat[-1,-1]
# 
# class(KDmat) <- "numeric"
# 
# #####################
# library(rentrez)
# entrez_db_summary("gene")
# x <- entrez_search(db="gene", term="4946")
# x$ids
# entrez_db_searchable("gene")
# sum <- entrez_summary(db = "gene", id = "105372840")
# sum$summary
# sum
# 
# 
# sym2eg <- org.Hs.egSYMBOL2EG
# str(sym2eg)


## ------------------------------------------------------------------ ##
##                          KD Overview                               ##
## ------------------------------------------------------------------ ##

library(tools)

# read KD data
cat(file = stderr(), "Processing KD data...\n")
KDdata <-
  read.csv(
    "./data/KDmat.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
pValP <-
  read.csv(
    "./data/pProx.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
pValD <-
  read.csv(
    "./data/pDist.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )

KDmat <- as.matrix(KDdata)
pValPmat <- as.matrix(pValP)
pValDmat <- as.matrix(pValD)

# save conditions as colnames and genes as rownames
colnames(KDmat) <- KDmat[1,]
rownames(KDmat) <- KDmat[, 1]
colnames(pValPmat) <- pValPmat[1,]
rownames(pValPmat) <- pValPmat[, 1]
colnames(pValDmat) <- pValDmat[1,]
rownames(pValDmat) <- pValDmat[, 1]

KDmat <- KDmat[-1,-1]
pValPmat <- pValPmat[-1,-1]
pValDmat <- pValDmat[-1,-1]

class(KDmat) <- "numeric"
class(pValPmat) <- "numeric"
class(pValDmat) <- "numeric"

# create overview matrix
# KDOverview <- KDmat
# KDOverview[KDOverview > 0] <- 1
# KDOverview[KDOverview < 0] <- 1
# KDOverview[is.na(KDOverview)] <- 0

# create matrix with 0s instead of NAs
KDmat0 <- KDmat
KDmat0[is.na(KDmat0)] <- 0
KDmat0 <- data.matrix(KDmat0)

# list with all conditions/KDs
conditions <- sort(colnames(KDmat))

# list with genes
genenames <- rownames(KDmat)

# create distance matrix
dist <- as.matrix(dist(abs(KDmat0)))
dist[dist == 0] <- NA

## ------------------------------------------------------------------ ##
##                          BigWigs                                   ##
## ------------------------------------------------------------------ ##

cat(file = stderr(), "Processing BigWig data...\n")

# get all bigwig files
bigwigs = list.files(path = "./data/IGV_vis_drops", pattern = ".bw")

# create matrix
bw_filemat <- matrix(, nrow = 0, ncol = 2)

# set colnames
colnames(bw_filemat) <- c("pos", "neg")

for (con in conditions) {
  x <- character(2)
  x <- NA
  
  bw_filemat <- rbind(bw_filemat, x)
  rownames(bw_filemat)[rownames(bw_filemat) == "x"] <- con
  
}

for (file in bigwigs) {
  x <- character(2)
  x <- NA
  
  filename = file_path_sans_ext(file)
  str = strsplit(filename, '_')
  # UTR2 / UTR3
  utr = str[[1]][1]
  # gene name
  gene = str[[1]][3]
  # strand (pos / neg)
  strand = str[[1]][4]
  
  # create row for each gene
  if (!(gene %in% rownames(bw_filemat))) {
    bw_filemat <- rbind(bw_filemat, x)
    rownames(bw_filemat)[rownames(bw_filemat) == "x"] <- gene
  }
  
  # add file for pos strand
  if (strand == "pos") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene, "pos"])) {
      
    }
    else
      bw_filemat[gene, "pos"] <- file
  }
  
  # add file for neg strand
  if (strand == "neg") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene, "neg"])) {
      
    }
    else
      bw_filemat[gene, "neg"] <- file
  }
}

# sort matrix by gene names
bw_filemat <- bw_filemat[order(rownames(bw_filemat)), ]


## ------------------------------------------------------------------ ##
##                          Gviz                                      ##
## ------------------------------------------------------------------ ##

library(Gviz)
library(org.Hs.eg.db)

#hg38db <- makeTxDbFromUCSC(genome = "hg38", tablename ="knownGene")
#saveDb(hg38db, file="hg38db_knownGene.sqlite")

cat(file = stderr(), "Loading TxDB...\n")

# load annotation db
txdb <- loadDb("./data/hg38db_refGene.sqlite")

# extract all genes as GRanges object
allgenes <- genes(txdb)

# get the entrez gene identifiers that are mapped to a gene symbol & save as list
sym2eg <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(sym2eg)
sym2eg_list <- as.list(sym2eg[mapped_genes])

# get 3' UTR data
utrs <- threeUTRsByTranscript(txdb)

# create data frame with gene & transcript ids
txBygene <- transcriptsBy(txdb, "gene")
gID <- rep(names(txBygene), elementNROWS(txBygene))
IDs <-
  data.frame(geneID = gID, txID = values(unlist(txBygene))[["tx_id"]])

# create genome axis track
gtrack <- GenomeAxisTrack()

# miRNA annotation
miRNA <- "sno_miRNA_hg38.bed"

# poly A site marker
polyApos <- "adj_contigs_pos_hg38.bed"
polyAneg <- "adj_contigs_neg_hg38.bed"


createLinkGO <- function (val) 
{
  sprintf("<a href=\"http://amigo.geneontology.org/amigo/term/%s\" target=\"_blank\" class=\"btn btn-primary\">%s</a>", 
          val, val)
}
