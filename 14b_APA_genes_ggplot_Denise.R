#! /usr/bin/r
#Script to get the plots of Fold Change of prox and dist contig for the entire gene
#two most signi contigs represent the whole gene behaviour
#based on individual list of genes for each KD

#on the plot OX = dist fold, OY = prox fold

library("gtools")
library("ggplot2")
library("grid")
library("RColorBrewer")

#List of files and sort in numerical "human" order
signi_files <- list.files(path="./data/14_APA_genes_whole_plus_dropouts")
signi_files <- mixedsort(signi_files)

#color palette
my_pal <- brewer.pal(11,"Spectral")[c(2,10)]

for (input in signi_files) {
  
  #extract KD name exclusively for column names in outputs
  kd_name <- as.character(data.frame(strsplit(input,"[_.]"))[2,])
  kd_name_to_plot <- as.character(data.frame(strsplit(input,"[.]"))[1,])
  
  #import file
  filepath <- file.path("./data/14_APA_genes_whole_plus_dropouts",input)
  import_data <- read.table(filepath, header = TRUE)
  
  #make a new table with relevanti info
  scatter <- cbind(import_data[,c(1,6,7)],APA_direction=as.character("Low_FC"), sort=4, Dir_idx=import_data[,8], stringsAsFactors=FALSE)
  
  #assign "Dist" for all distalisations (Dir_idx < 0)
  scatter[which(scatter$Dir_idx <= - 0.5),4] <- "Dist"
  #assign "Prox" for proximisation ((Dir_idx > 0))
  scatter[which(scatter$Dir_idx >= 0.5),4] <- "Prox"
  
  #assign sort "1" for distalizations
  scatter[which(scatter$APA_direction == "Dist"),5] <- 1
  #assign sort "2" for proximizations
  scatter[which(scatter$APA_direction == "Prox"),5] <- 2
  #assign sort "3" for Low_FC
  scatter[which(scatter$APA_direction == "Low_FC"),5] <- 3
  
  #order by "sort"
  scatter <- scatter[order(scatter$sort, decreasing=TRUE),]
  
  #make an annotation for ndist and nprox ("grid" package)
  dist <- grobTree(textGrob(paste("dist",sum(scatter[,4]=="Dist")), x=0.65,  y=0.05, hjust=0,
                            gp=gpar(col=my_pal[1], fontsize=36)))
  prox <- grobTree(textGrob(paste("prox",sum(scatter[,4]=="Prox")), x=0.05,  y=0.95, hjust=0,
                            gp=gpar(col=my_pal[2], fontsize=36)))
  
  
  #make grid with OX and OY
  ox <- grobTree(linesGrob(c(0,1), c(0,1), gp=gpar(lwd=1.5)))
  
  color <- scatter$APA_direction
  scatter <- scatter[,2:3]
  #make a file name for each plot
  outpath <- file.path("./data/scatterplots", paste(kd_name_to_plot, ".jpg", sep=""))
  
  #make a plot for each KD (manually), for significant contigs only
  #where a fold=f(delta_length)
  
  ggplot(scatter, aes(Fold_dist,Fold_prox,colour=color), na.rm=TRUE)+
    annotation_custom(ox)+  #adds line
    geom_point(size=3.4)+ #size of the dots
    scale_color_manual(values=c(my_pal[1],"black", my_pal[2]),name="APA direction")+  #color of the dots and legend name
    xlim(-6,6)+
    ylim(-10,10)+
    labs(x="Log2(FC), distal PAS", y="Log2(FC), proximal PAS")+  #labels
    ggtitle(kd_name)+ #title
    #theme adjustments
    theme(plot.title = element_text(size=18, face="bold", vjust=2), 
          panel.background = element_rect(fill = 'gray98'),
          axis.text.x=element_text(size=30, vjust=0.5),
          axis.text.y=element_text(size=30, vjust=0.5),
          axis.title.x = element_text(size=30, vjust=-0.35),
          axis.title.y = element_text(size=30, vjust=0.35))+
    annotation_custom(dist)+
    annotation_custom(prox)
  
  
  ggsave(filename=outpath,width=20,height=18,units="cm")
  
}








