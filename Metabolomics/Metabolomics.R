library("metabolomics")
library("vioplot")
library("ggplot2")
setwd("/Users/elvirakinzina/src/Bioinfo/Metabolomics")

read.file <- function(file) {
  raw_data <- read.csv(paste("GC/final/",toString(file),sep = ""), sep = '\t', na.strings = c("NA", "NaN"))
  if (file!="GC Liver final_No IS") {
     raw_data$putative.Compound <- gsub(";","",raw_data$putative.Compound)
  }
  data_no_missing <- raw_data[which(rowSums(is.na(raw_data))<dim(raw_data)[2]-1),]
  return(data_no_missing)
}

overall.boxplot <- function(file, data) {
  png(file=file, width = 1000, height = 400)
  new_data <- data[2:dim(data)[2]]
  column.names <- gsub("_","",substr(colnames(new_data),1,7))
  boxplot(new_data, names = c(column.names), ylim = range(0,200000))
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

plot.sign.comp.boxplot <- function(compound.data,j) {
  print(compound.data$putative.Compound[j])
  filename = paste("No permutation/Boxplots/Sign met/", toString(compound.data$putative.Compound[j]), ".png", sep="")
  png(file=filename)
  file <- read.table("Files.txt", sep = "\n")
  
  d <- data.frame(x = numeric(0), col = numeric(0), mean = numeric(0), se = numeric(0))
  tissues <- NULL
  for (i in 1:dim(file)[1]) {
    data <- read.file(file$V1[i])
    compound <- data$putative.Compound[j]
    data <- data[,2:dim(data)[2]]
    has.na <- apply(is.na(data),1,sum)>0
    data <- data[!has.na,]
    data <- data[, !colnames(data) %in% c("Wt_91_1._cer")]
    dh <- data[grep("hum",colnames(data))]
    dw <- data[grep("Wt",colnames(data))]
    dH <- as.numeric(as.matrix(dh[j,]))
    dW <- as.numeric(as.matrix(dw[j,]))
    mean.dh <- mean(dH)
    mean.dw <- mean(dW)
    se.h <- sd(dH)/sqrt(length(dH))
    se.w <- sd(dW)/sqrt(length(dW))
    d <- rbind(d, data.frame(x = i, col = 1, mean = mean.dh, se = se.h))
    d <- rbind(d, data.frame(x = i, col = 2, mean = mean.dw, se = se.w))
    tissues[i] <- strsplit(toString(file$V1[i]),"_")[[1]][2]
  }
  colors <- c("red","blue")[d$col]
ggplot(d, aes(x=x, y=mean)) + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, color = colors) +
    geom_point(color=colors) + labs(list(title = compound, x = "type of tissue", y = "mean concentration")) +
    scale_x_continuous(breaks = 1:9, labels = tissues)
  ggsave(filename, device = "png")
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
    sorted.indexes <- sort(as.vector(content[!is.na(content)]), decreasing = FALSE, index.return=TRUE, na.last = NA)$ix
    for (k in 1:length(sorted.indexes)) {
      write(paste(toString(data$putative.Compound[index[sorted.indexes[k]]]),
                  toString(sub("\\.",",",round(content[sorted.indexes[k]], digits = 3))),sep=";"), 
            file = file_print, append=TRUE)
      #print(data$putative.Compound[index[sorted.indexes[k]]])
    }
  }
}

perm <- function(data) {
  shuffled.data <- sample(data[,2:dim(data)[2]],dim(data)[2]-1)
  human <- shuffled.data[,1:dim(data[grep("hum",colnames(data))])[2]]
  wild <- shuffled.data[,(dim(data[grep("hum",colnames(data))])[2]+1):dim(shuffled.data)[2]]
  N <- min(dim(human)[1], dim(wild)[1])
  p.value <- NULL
  for (j in 1:N) {
    result <- t.test(as.numeric(as.matrix(human[j,])), as.numeric(as.matrix(wild[j,])))
    p.value[j] <- result$p.value
  }
  n.sign.met <- sum(p.value<0.05)
  return(n.sign.met)
}

main <- function() {
  files <- read.table("Files.txt", sep = "\n")
  
  for (i in 1:dim(files)[1]) {
    file <- files$V1[i]
    print(file)
    data <- read.file(file)
    #data <- data[which(grepl("FAME",data$putative.Compound)==FALSE),]
    Metabolites <- data$putative.Compound
    #pca.result.new <- prcomp(as.matrix(data[2:dim(data)[2]]), na.action=na.omit(data), scale = FALSE)
    pca.result <- pca(data)
    summary(pca.result)
    filename = paste("No permutation/PCA/", toString(file), ".png", sep="")
    png(file=filename, width = 1000, height = 600)
    slplot(pca.result)
    dev.off()
    print(dim(data)[2])
    #remove outliers
    if (i==1) {
      data <- data[, !colnames(data) %in% c("Wt_91_1._cer")]
    }
    if (i==2) {
      data <- data[, !colnames(data) %in% c("Wt_135_2._cor")]
    }
    if (i==3) {
      data <- data[, !colnames(data) %in% c("hum_120_6._kid")]
      #Wt_91_6._kid
    }
    if (i==4) {
      data <- data[, !colnames(data) %in% c("hum_105_3._hea.1")]
      #hum_98_3._hea
    }
    if (i==5) {
      data <- data[, !colnames(data) %in% c("hum_120_9._mus.1")]
      ##hum_121_9._mus.1
    }
    if (i==7) {
      data <- data[, !colnames(data) %in% c("hum_101_7._spl")]
      #hum_120_7._spl
      ##hum_153_7_spl
    }
    if (i==8) {
      data <- data[, !colnames(data) %in% c("hum_101_8._tes")]
      #hum_120_8._tes
    }
    if (i==9) {
      data <- data[, !colnames(data) %in% c("Wt_91_5._liv")]
      #hum_101_5._liv
    }
    print("new")
    print(dim(data)[2])
    data_hum <- data[grep("hum",colnames(data))]
    data_hum <- data_hum[which(rowSums(is.na(data_hum))<dim(data_hum)[2]),]
    data_wt <- data[grep("Wt",colnames(data))]
    data_wt <- data_wt[which(rowSums(is.na(data_wt))<dim(data_wt)[2]),]
    
    file_print = paste("No permutation/", "/Boxplots/", toString(file)," boxplot.png",sep="")
      #overall.boxplot(file_print, data)
      
    adjusted.p.value <- NULL
    just.p.value <- NULL
    sign.p.value <- NULL
    just.sign.p.value <- NULL
    indexes <- NULL
    just.indexes <- NULL
    
    N <- min(dim(data_hum)[1], dim(data_wt)[1])
    for (j in 1:N) {
      result <- t.test(as.numeric(as.matrix(data_hum[j,])), as.numeric(as.matrix(data_wt[j,])))
      just.p.value[j] <- result$p.value
    }
    adjusted.p.value <- p.adjust(just.p.value, method="BH")
    d <- data.frame(metabolites = character(0), p.val = numeric(0), adj.p.val = numeric(0))
    for (l in 1:length(just.p.value)) {
      d <- rbind(d, data.frame(metabolites = Metabolites[order(just.p.value, decreasing = FALSE)][l], p.val = sort(just.p.value, decreasing = FALSE)[l], adj.p.val = sort(adjusted.p.value, decreasing = FALSE)[l]))
    }
    file_print = paste("No permutation/P-value/", toString(file)," compounds p-value.csv",sep="")
    write.table(d, file = file_print, sep = ";")
    n.sign.met <- sum(just.p.value<0.05)
    
    #p-value distribution
    filename = paste("No permutation/P-value distribution/", toString(files$V1[i]), "p-value distr.png", sep="")
    png(file=filename)
    hist(as.numeric(just.p.value), main = "P-value distribution", xaxt='n', breaks = length(as.numeric(just.p.value)[!is.na(just.p.value)])/7)
    axis(side=1, at=seq(0,1,0.05), labels=seq(0,1,0.05))
    dev.off()
    
    #permutations
    # permutation.number <- 1000
    # number.of.sign.met <- NULL
    # for (m in 1:permutation.number) {
    #   print(m)
    #   number.of.sign.met[m] <- perm(data)
    # }
    # write.csv(number.of.sign.met, file = paste("Permutation/", toString(file), sep=""))
    # filename = paste("Permutation/", toString(files$V1[i]), ".png", sep="")
    # png(file=filename)
    # hist(as.numeric(number.of.sign.met), main = paste("Number of significant metabolites for ", toString(permutation.number), " permutations \n (", toString(n.sign.met), " without permutations)", sep = ""), breaks = length(as.numeric(number.of.sign.met)[!is.na(number.of.sign.met)])/50)
    # dev.off()

    #file.p.value = paste("No permutation/P-value distribution/", toString(file), " p-value.png", sep="")
    #hist.p.value.to.file(just.p.value, N, file.p.value, "Distribution of p-value", "p-value")
  }
}
main()
