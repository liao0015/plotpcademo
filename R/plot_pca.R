#' plot PCA
#'
#' This function loads files and plot PCA as output 
#' @param input1 experimental design input
#' @param proteinGroups protein groups input
#' @export

plot_pca<-function(input1, proteinGroups){


# read in experimental desgin
#ExperimentalDesign<-read.delim("inst/extdata/Experimental_Design.txt",header=TRUE)
#sort according to the file names
ExperimentalDesign<-input1[order(input1$Old.Sample.Name),]

# read in full proteinGroups.txt file
#proteinGroups<-read.delim("inst/extdata/proteinGroups.txt",row.names=1,header=TRUE)

# read in configuration from html in the future
  
  proteinGroups_Reversed<-proteinGroups[proteinGroups$Reverse=="+",]
  proteinGroups_Contaminant<-proteinGroups[proteinGroups$Contaminant=="+",]
  proteinGroups_Only.identified.by.site<-proteinGroups[proteinGroups$Only.identified.by.site=="+",]
  
  # count how many for the three reverse/contaminat/identified.by.site  
  proteinGroups_count<-nrow(proteinGroups)
  proteinGroups_count_Reversed<-nrow(proteinGroups_Reversed)
  proteinGroups_count_Contaminant<-nrow(proteinGroups_Contaminant)  
  proteinGroups_count_Only.identified.by.site<-nrow(proteinGroups_Only.identified.by.site)
  
  # remove all the rows marked as reverse/contaminat/identified.by.site
  proteinGroups_filtered<-proteinGroups[proteinGroups$Reverse!="+" & proteinGroups$Contaminant!="+" & proteinGroups$Only.identified.by.site!="+",]
  proteinGroups_filtered_count<-nrow(proteinGroups_filtered)
  
  proteinGroups_filtering<-(list(filtered=proteinGroups_filtered,
              proteinGroups_count=proteinGroups_count,
              proteinGroups_count_Reversed=proteinGroups_count_Reversed,
              proteinGroups_count_Contaminant=proteinGroups_count_Contaminant,
              proteinGroups_count_Only.identified.by.site=proteinGroups_count_Only.identified.by.site,
              proteinGroups_filtered_count=proteinGroups_filtered_count    
  ))
  
# Step 1: do filtering, to remove rows with contaminant/revers/identified by id
proteinGroups_filtered<-proteinGroups_filtering$filtered

#The filtering summary can also be reported out to the user end
#proteingroups before filtering:
#proteinGroups_filtering$proteinGroups_count
#proteingroups after filtering:
#proteinGroups_filtering$proteinGroups_filtered_count
#proteingroups marked as reversed
#proteinGroups_filtering$proteinGroups_count_Reversed
#proteingroups  marked as contaminant
#proteinGroups_filtering$proteinGroups_count_Contaminant
#proteingroups marked as Only identified by site
#proteinGroups_filtering$proteinGroups_count_Only.identified.by.site


# keep the LFQ columns as example for downstream analysis
proteinGroups_headers<-colnames(proteinGroups_filtered)
proteinGroups_filtered_LFQ_intensity<-proteinGroups_filtered[,grep("LFQ.Intensity.",proteinGroups_headers)]

# simplify LFQ column names
colnames(proteinGroups_filtered_LFQ_intensity)<-gsub("LFQ.Intensity.", "", colnames(proteinGroups_filtered_LFQ_intensity))
# re-arrange the columns, according to the names
proteinGroups_filtered_LFQ_intensity<-proteinGroups_filtered_LFQ_intensity[,order(colnames(proteinGroups_filtered_LFQ_intensity))]



# Step 2: very basic data process: normalization

# replace the 0 with NaN
proteinGroups_filtered_LFQ_intensity[proteinGroups_filtered_LFQ_intensity==0]<-NaN
# take log10
proteinGroups_filtered_LFQ_intensity_log10<-log10(proteinGroups_filtered_LFQ_intensity)

# filter out the entrieséproteins which have larger number of NA counting
# the threshold cold be easily setup by QXX starndard, by using ncol(matrix)*Q
# here we start from with 0 tolerannce of NA (Q100) as a test
proteinGroups_filtered_LFQ_intensity_log10_Q100<-proteinGroups_filtered_LFQ_intensity_log10[which(apply(proteinGroups_filtered_LFQ_intensity_log10,1,function(x)(sum(is.na(x)))<1)),]

# simple displace the density sum of each sample, to see if there is any sample effect, Q100 proteins are very high abundant ones, which cann tell 
boxplot(proteinGroups_filtered_LFQ_intensity_log10_Q100,col = "blue", xlab="Samples",cex.axis=0.5, las=2,main="LFQ Intensity Profile")

# do scaling, keep in mind that the scale function in R is scaling by column
proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled<-scale(proteinGroups_filtered_LFQ_intensity_log10_Q100)

# double check the effect of scaling by comparing the box-plot with the above
boxplot(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled,col = "blue", xlab="Samples",cex.axis=0.5, las=2,main="LFQ Intensity Profile Normalzed Across Samples")

#further check the scaling effect:
plot(apply(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled, 1, mean), ylim=c(-5,5),main="Mean of Proteins", xlab="Proteins Index", ylab = "Average Protein Abundance of proteins Across Samples")
barplot(apply(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled, 2, mean), ylim=c(-5,5),main="Mean of Samples", xlab= "Samples", ylab="Average Normalized LFQ Intensenties")

# scaling of each protein
proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled<-t(scale(t(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled)))
plot(apply(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled, 1, mean), ylim=c(-5,5),main="Mean of Proteins",xlab="Proteins Index", ylab = "Average Protein Abundance of proteins Across Samples")
barplot(apply(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled, 2, mean), ylim=c(-5,5),main="Mean of Samples", xlab= "Samples", ylab="Average Normalized LFQ Intensenties")

# proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled will be the final processed data to be used

# however, since some funtions have built in "scaling", you can choose whcihever staring point


# Step 3: PCA plot, with screeplot 
pca.Inensity <- prcomp(t(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled), scale.=TRUE, center = TRUE) 
plot(pca.Inensity$x,col = ExperimentalDesign$Groups, pch=15, main="PCA plot")
calibrate::textxy(pca.Inensity$x[,1],pca.Inensity$x[,2], labs = rownames(pca.Inensity$x))

  pca.Inensity <- prcomp(t(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled), scale.=TRUE, center = TRUE) 
  sd <- pca.Inensity$sdev
  loadings <- pca.Inensity$rotation
  rownames(loadings) <- rownames(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled)
  scores <- pca.Inensity$x
  var <- sd^2
  var.percent <- var/sum(var) * 100
  barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray", main="Percent of Variance")
  abline(h=1/nrow(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled)*100, col="red")


}







