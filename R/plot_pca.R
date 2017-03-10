#' plot PCA
#'
#' This function loads files and plot PCA as output 
#' @export

plot_p<-function(){

# two files as input:

# read in experimental desgin
ExperimentalDesign<-read.delim("inst/extdata/Experimental_Design.txt",header=TRUE)
#sort according to the file names
ExperimentalDesign<-ExperimentalDesign[order(ExperimentalDesign$Old.Sample.Name),]

# read in full proteinGroups.txt file
proteinGroups<-read.delim("inst/extdata/proteinGroups.txt",row.names=1,header=TRUE)

# read in configuration from html in the future



# Step 1: do filtering, to remove rows with contaminant/revers/identified by id
proteinGroups_filtering<-proteinGroups_filter(proteinGroups)
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
# pca.Inensity <- prcomp(t(proteinGroups_filtered_LFQ_intensity_log10_Q100_scaled.scaled), scale.=TRUE, center = TRUE) 
# plot(pca.Inensity$x,col = ExperimentalDesign$Groups, pch=15, main="PCA plot")
# calibrate::textxy(pca.Inensity$x[,1],pca.Inensity$x[,2], labs = rownames(pca.Inensity$x))

}









