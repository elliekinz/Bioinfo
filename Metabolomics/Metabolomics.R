library("metabolomics")
library("vioplot")
setwd("/Users/elvirakinzina/src/Bioinfo/Metabolomics")
files <- read.table("Files.txt", sep = "\n")
for (i in 1:dim(files)[1]) {
  i <- 9
  data <- read.csv(paste("GC/",toString(files$V1[i]),".txt",sep = ""), sep = '\t',na.strings = "NaN")
  data_no_missing <- data[which(rowSums(is.na(data))<dim(data)[2]-1),]
  if (files$V1[i]!="GC Liver final_No IS") {
    data_no_missing$putative.Compound <- gsub(";","",data_no_missing$putative.Compound)
  }
  
  #data_no_na <- data_no_missing[which(rowSums(is.na(data_no_missing[grep("hum",colnames(data_no_missing))][j,]))<dim(data_no_missing[grep("hum",colnames(data_no_missing))][j,])[2]-1),]
  #data_no_na <- data_no_na[which(rowSums(is.na(data_no_na[grep("Wt",colnames(data_no_na))][j,]))<dim(data_no_na[grep("Wt",colnames(data_no_na))][j,])[2]-1),]
  adjusted.p.value <- NULL
  just.p.value <- NULL
  sign.p.value <- NULL
  just.sign.p.value <- NULL
  indexes <- NULL
  just.indexes <- NULL
  print(files$V1[i])
  for (j in 2:dim(data)[1]) {
    if (rowSums(is.na(data_no_missing[grep("hum",colnames(data_no_missing))][j,])) < 
        dim(data_no_missing[grep("hum",colnames(data_no_missing))][j,])[2]-1) {
      # filename = paste("Distribution/", toString(files$V1[i]),"/Boxplots/box", j, ".png", sep="")
      # png(file=filename)
      # boxplot(as.numeric(data_no_missing[grep("hum",colnames(data_no_missing))][j,]),
      #         as.numeric(data_no_missing[grep("Wt",colnames(data_no_missing))][j,]), 
      #         names = c(paste("hum",substr(toString(data_no_missing[j,]$putative.Compound),1,10)),
      #                   paste("Wt",substr(toString(data_no_missing[j,]$putative.Compound),1,10))))
      # dev.off() 
      result <- t.test(as.numeric(data_no_missing[grep("hum",colnames(data_no_missing))][j,]),
                       as.numeric(data_no_missing[grep("Wt",colnames(data_no_missing))][j,]))
      just.p.value[j] <- result$p.value
      adjusted.p.value[j] <- p.adjust(just.p.value[j], method="bonferroni")
      
      if ((just.p.value[j] < 5e-2) & (!is.null(data_no_missing$putative.Compound[j]))) {
        if (!is.na(data_no_missing$putative.Compound[j])) {
          just.sign.p.value[length(just.sign.p.value)+1] <- just.p.value[j]
          just.indexes[length(just.indexes)+1] <- j
        }
      }
      
      if ((adjusted.p.value[j] < 5e-2) & (!is.null(data_no_missing$putative.Compound[j]))) {
        if (!is.na(data_no_missing$putative.Compound[j])) {
          sign.p.value[length(sign.p.value)+1] <- adjusted.p.value[j]
          indexes[length(indexes)+1] <- j
        }
      }
    }
  }
  
  if (!is.null(sign.p.value)) {
    file_print = paste("Print/",toString(files$V1[i]),".csv",sep="")
    sorted.indexes <- sort(as.vector(sign.p.value[!is.na(sign.p.value)]), decreasing = TRUE, index.return=TRUE, na.last = NA)$ix
    for (k in 1:length(sorted.indexes)) {
      write(paste(toString(data_no_missing$putative.Compound[indexes[sorted.indexes[k]]]),
                  toString(sub("\\.",",",round(sign.p.value[sorted.indexes[k]], digits = 3))),sep=";"), 
            file = file_print, append=TRUE)
    }
  }
  
  if (!is.null(just.sign.p.value)) {
    file_print = paste("Print/",toString(files$V1[i])," p-value.csv",sep="")
    just.sorted.indexes <- sort(as.vector(just.sign.p.value[!is.na(just.sign.p.value)]), decreasing = TRUE, index.return=TRUE, na.last = NA)$ix
    for (k in 1:length(just.sorted.indexes)) {
      write(paste(toString(data_no_missing$putative.Compound[just.indexes[just.sorted.indexes[k]]]),
                  toString(sub("\\.",",",round(just.sign.p.value[just.sorted.indexes[k]], digits = 3))),sep=";"), 
            file = file_print, append=TRUE)
    }
  }
  
  file.p.value = paste("Print/P_value/", toString(files$V1[i]), " adjusted p-value.png", sep="")
  png(file=file.p.value)
  hist(as.numeric(adjusted.p.value), main = 'Distribution adjusted of p-value', xlab = 'adjusted p-value',
       breaks = length(as.numeric(adjusted.p.value)[!is.na(adjusted.p.value)])/2, xaxt='n')
  axis(side=1, at=seq(0,1,0.05), labels=seq(0,1,0.05))
  dev.off()

  file.p.value = paste("Print/P_value/", toString(files$V1[i]), " p-value.png", sep="")
  png(file=file.p.value)
  hist(as.numeric(just.p.value), main = 'Distribution of p-value', xlab = ' p-value',
       breaks = length(as.numeric(just.p.value)[!is.na(just.p.value)])/2, xaxt='n')
  axis(side=1, at=seq(0,1,0.05), labels=seq(0,1,0.05))
  dev.off()
  
}
