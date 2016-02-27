library("metabolomics")
library("vioplot")
setwd("/Users/elvirakinzina/src/Bioinfo/Metabolomics")

read.file <- function(file) {
  raw_data <- read.csv(paste("GC/",toString(file),".txt",sep = ""), sep = '\t',na.strings = "NaN")
  data_no_missing <- raw_data[which(rowSums(is.na(raw_data))<dim(raw_data)[2]-1),]
  if (file!="GC Liver final_No IS") {
    data_no_missing$putative.Compound <- gsub(";","",data_no_missing$putative.Compound)
  }
  return(data_no_missing)
}

overall.boxplot <- function(file, data) {
  png(file=file, width = 1400, height = 700)
  new_data <- data[2:dim(data)[2]]
  column.names <- gsub("_","",substr(colnames(new_data),1,7))
  boxplot(new_data, names = c(column.names))
  dev.off()
}

plot.boxplot <- function(data) {
  filename = paste("Distribution/", toString(files$V1[i]),"/Boxplots/box", j, ".png", sep="")
  png(file=filename)
  boxplot(as.numeric(data[grep("hum",colnames(data))][j,]),
          as.numeric(data[grep("Wt",colnames(data))][j,]),
          names = c(paste("hum",substr(toString(data[j,]$putative.Compound),1,10)),
                    paste("Wt",substr(toString(data[j,]$putative.Compound),1,10))))
  dev.off()
}

hist.p.value.to.file <- function(content, n, file, title, axis.name) {
  png(file=file, width = 1000, height = 600)
  hist(as.numeric(content), main = title, xlab = axis.name, xaxt='n', breaks = length(as.numeric(content)[!is.na(content)])/2)
  axis(side=1, at=seq(0,1,0.05), labels=seq(0,1,0.05))
  dev.off()
}

write.sign.comp.to.file <- function(index, content, data, file_print) {
  if (!is.null(content)) {
    sorted.indexes <- sort(as.vector(content[!is.na(content)]), decreasing = TRUE, index.return=TRUE, na.last = NA)$ix
    for (k in 1:length(sorted.indexes)) {
      write(paste(toString(data$putative.Compound[index[sorted.indexes[k]]]),
                  toString(sub("\\.",",",round(content[sorted.indexes[k]], digits = 3))),sep=";"), 
            file = file_print, append=TRUE)
    }
  }
}

main <- function(permutation) {
  files <- read.table("Files.txt", sep = "\n")
  
  for (i in 1:dim(files)[1]) {
    file <- files$V1[i]
    data <- read.file(file)
    if (permutation != 0) {
      file_print = paste("Permutation", toString(permutation), "/Boxplots/", toString(file)," boxplot.png",sep="")
    }
    else {
      file_print = paste("No permutation/", "/Boxplots/", toString(file)," boxplot.png",sep="")
    }
    overall.boxplot(file_print, data)
    N <- dim(data)[1]
    #data_no_na <- data_no_missing[which(rowSums(is.na(data_no_missing[grep("hum",colnames(data_no_missing))][j,]))<dim(data_no_missing[grep("hum",colnames(data_no_missing))][j,])[2]-1),]
    #data_no_na <- data_no_na[which(rowSums(is.na(data_no_na[grep("Wt",colnames(data_no_na))][j,]))<dim(data_no_na[grep("Wt",colnames(data_no_na))][j,])[2]-1),]
    adjusted.p.value <- NULL
    just.p.value <- NULL
    sign.p.value <- NULL
    just.sign.p.value <- NULL
    indexes <- NULL
    just.indexes <- NULL
    shuffed.data <- sample(data[,2:dim(data)[2]],dim(data)[2]-1)
    for (j in 2:dim(data)[1]) {
      if (rowSums(is.na(data[grep("hum",colnames(data))][j,])) < 
          dim(data[grep("hum",colnames(data))][j,])[2]-1) {
        # plot.boxplot(data)
        if (permutation != 0) {
          result <- t.test(as.numeric(shuffed.data[,1:dim(data[grep("hum",colnames(data))])[2]][j,]),
                           as.numeric(shuffed.data[,(dim(data[grep("hum",colnames(data))])[2]+1):dim(shuffed.data)[2]][j,]))
        }
        else {
          result <- t.test(as.numeric(data[grep("hum",colnames(data))][j,]),
                           as.numeric(data[grep("Wt",colnames(data))][j,]))
        }
        just.p.value[j] <- result$p.value
      }
      else {
        just.p.value[j] <- NA
      }
    }
    
    just.p.value <- just.p.value[!is.na(just.p.value)]
    adjusted.p.value <- p.adjust(just.p.value, method="bonferroni")
    
    threshold <- 5e-2
    for (k in 1:length(just.p.value)) {
      if ((just.p.value[k] < threshold) & (!is.null(data$putative.Compound[k]))) {
        if (!is.na(data$putative.Compound[k])) {
          just.sign.p.value[length(just.sign.p.value)+1] <- just.p.value[k]
          just.indexes[length(just.indexes)+1] <- k
        }
      }
      if ((adjusted.p.value[k] < threshold) & (!is.null(data$putative.Compound[k]))) {
        if (!is.na(data$putative.Compound[k])) {
          sign.p.value[length(sign.p.value)+1] <- adjusted.p.value[k]
          indexes[length(indexes)+1] <- k
        }
      }
    }
    
    if (permutation != 0) {
      file_print = paste("Permutation", toString(permutation), "/", toString(file)," compounds adjusted p-value.csv",sep="")
    }
    else {
      file_print = paste("No permutation/", toString(file)," compounds adjusted p-value.csv",sep="")
    }
    write.sign.comp.to.file(indexes, sign.p.value, data, file_print)
    
    if (permutation != 0) {
      file.p.value = paste("Permutation", toString(permutation), "/P-value distribution/", toString(file), " p-value.png", sep="")
    }
    else {
      file.p.value = paste("No permutation/P-value distribution/", toString(file), " p-value.png", sep="")
    }
    hist.p.value.to.file(just.p.value, N, file.p.value, "Distribution of p-value", "p-value")
  }
}

main(0)
main(1)
main(2)
main(3)