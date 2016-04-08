library(tools)

setwd("/Users/Denise/Documents/Master/IMBEI/data/14_APA_genes_whole_plus_dropouts")

temp = list.files(pattern=".txt")

# create empty matrix with a column for each condition
KD_matrix <- matrix(,nrow = 0, ncol = length(temp))

# create list with all conditions/KDs
all_conditions <- list()

for (file in temp) {
  filename = file_path_sans_ext(file)    
  str = strsplit(filename, '_')
  condition = str[[1]][2]
  all_conditions[[length(all_conditions)+1]] <- condition
}

# set colnames in matrix to conditions
colnames(KD_matrix) <- all_conditions

# read all condition files & create row in matrix for every gene, then fill in dir_idx for all affected genes
for (file in temp) {
  
  data <- read.delim(file, stringsAsFactors = TRUE)
  
  for (gene in data$Gene_symbol) {
    
    x <- character(length(temp))
    x <- NA
    
    if(!(gene %in% rownames(KD_matrix))) {
      KD_matrix <- rbind(KD_matrix,x)
      rownames(KD_matrix)[rownames(KD_matrix)=="x"] <- gene
    }
    
    filename = file_path_sans_ext(file)    
    str = strsplit(filename, '_')
    condition = str[[1]][2]
    
    #condition <- filename
    
    # extract dir index
    dir <- data[data$Gene_symbol == gene, "Dir_idx"]
    
    KD_matrix[gene, condition] <- dir
  }
  
}

# sort matrix by gene names (= rownames)
KD_matrix <- KD_matrix[order(rownames(KD_matrix)), ] 

# delete rows with all NAs (i.e. genes that are not affected in any condition)
KD_matrix_relevant <- KD_matrix[rowSums(is.na(KD_matrix)) != ncol(KD_matrix)-1,]

setwd("/Users/Denise/Documents/Master/IMBEI/data/")

# save matrix in file
write.csv2(KD_matrix_relevant, file = "KDmat.csv")

# list with genes
genenames <- rownames(KD_matrix_relevant)