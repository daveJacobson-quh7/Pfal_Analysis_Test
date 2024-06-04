library(readxl)

time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()

timestamp <- as.Date(date_now, format = "%Y-%m-%d")



args <- commandArgs(trailingOnly = TRUE)


gatk_country_predictionDF <- read.table(args[1], sep = "\t", header =T)
gatk_region_predictionDF <- read.table(args[2], sep = "\t", header =T)

bcftools_country_predictionDF <- read.table(args[3], sep = "\t", header =T)
bcftools_region_predictionDF <- read.table(args[4], sep = "\t", header =T)

gatk_barcodeCSV <- read.table(args[5], sep =",", header = T)
bcftools_barcodeCSV <- read.table(args[6], sep =",", header = T)
metadataDF <- read_excel(args[7])

gatk_continent_predictionDF <- read.table(args[8], sep = "\t", header =T)
bcftools_continent_predictionDF <- read.table(args[9], sep = "\t", header =T)

# gatk_country_predictionDF <- gatk_country_predictionDF[,c(1:10)]
# gatk_region_predictionDF <- gatk_region_predictionDF[,c(1:10)]

# bcftools_country_predictionDF <- bcftools_country_predictionDF[,c(1:10)]
# bcftools_region_predictionDF <- bcftools_region_predictionDF[,c(1:10)]

gatk_country_predictionDF$Country_Prediction_1<- paste(gatk_country_predictionDF$Prediction_1, " (", round(gatk_country_predictionDF$Prob_1, digits = 5),")", sep = "")
gatk_country_predictionDF$Country_Prediction_2 <- paste(gatk_country_predictionDF$Prediction_2, " (", round(gatk_country_predictionDF$Prob_2, digits = 5),")", sep = "")
gatk_country_predictionDF$Country_Prediction_3 <- paste(gatk_country_predictionDF$Prediction_3, " (", round(gatk_country_predictionDF$Prob_3, digits = 5),")", sep = "")

gatk_region_predictionDF$Region_Prediction_1<- paste(gatk_region_predictionDF$Prediction_1, " (", round(gatk_region_predictionDF$Prob_1, digits = 5),")", sep = "")
gatk_region_predictionDF$Region_Prediction_2 <- paste(gatk_region_predictionDF$Prediction_2, " (", round(gatk_region_predictionDF$Prob_2, digits = 5),")", sep = "")
gatk_region_predictionDF$Region_Prediction_3 <- paste(gatk_region_predictionDF$Prediction_3, " (", round(gatk_region_predictionDF$Prob_3, digits = 5),")", sep = "")

gatk_continent_predictionDF$continent_Prediction_1<- paste(gatk_continent_predictionDF$Prediction_1, " (", round(gatk_continent_predictionDF$Prob_1, digits = 5),")", sep = "")
gatk_continent_predictionDF$continent_Prediction_2 <- paste(gatk_continent_predictionDF$Prediction_2, " (", round(gatk_continent_predictionDF$Prob_2, digits = 5),")", sep = "")
gatk_continent_predictionDF$continent_Prediction_3 <- paste(gatk_continent_predictionDF$Prediction_3, " (", round(gatk_continent_predictionDF$Prob_3, digits = 5),")", sep = "")



sampleOrder <- gatk_barcodeCSV$Sample
newVec <- c()
# for(i in sampleOrder){
#     newVar <- strsplit(i, split = "_")[[1]][1]
#     print(newVar)
#     newVec <- c(newVec, newVar)
#     print(newVec)
# }
#WGS samples have underscores in name, so have to adjust code
for(i in sampleOrder){
    newVar <- strsplit(i, split = "_sWGA")[[1]][1]
    print(newVar)
    newVec <- c(newVec, newVar)
    print(newVec)
}
gatk_country_predictionDF$Sample <- newVec
# gatk_country_predictionDF$Sample <- newVarL001","",gatk_country_predictionDF$Sample)
print(gatk_country_predictionDF$Sample)
gatk_region_predictionDF$Sample <- newVec
gatk_continent_predictionDF$Sample <- newVec



#Sample | pred1 |pred2 | pred3 | hets | missing | masked
print("get colnames")
# print(colnames(gatk_country_predictionDF))

gatk_country_predictionDF <- gatk_country_predictionDF[,c(15,12,13,14,8,9,10)]
gatk_region_predictionDF <- gatk_region_predictionDF[,c(15,12,13,14,8,9,10)]
gatk_continent_predictionDF <- gatk_continent_predictionDF[,c(15,12,13,14,8,9,10)]
# print(gatk_country_predictionDF)
print(("myMeta"))
print(head(metadataDF))
firstMerge <- merge(gatk_country_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)
secondMerge <- merge(gatk_region_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)
fifthMerge <- merge(gatk_continent_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)

print("afterMerge")
print(head(secondMerge))

bcftools_country_predictionDF$Country_Prediction_1<- paste(bcftools_country_predictionDF$Prediction_1, " (", round(bcftools_country_predictionDF$Prob_1, digits = 5),")", sep = "")
bcftools_country_predictionDF$Country_Prediction_2 <- paste(bcftools_country_predictionDF$Prediction_2, " (", round(bcftools_country_predictionDF$Prob_2, digits = 5),")", sep = "")
bcftools_country_predictionDF$Country_Prediction_3 <- paste(bcftools_country_predictionDF$Prediction_3, " (", round(bcftools_country_predictionDF$Prob_3, digits = 5),")", sep = "")

bcftools_region_predictionDF$Region_Prediction_1<- paste(bcftools_region_predictionDF$Prediction_1, " (", round(bcftools_region_predictionDF$Prob_1, digits = 5),")", sep = "")
bcftools_region_predictionDF$Region_Prediction_2 <- paste(bcftools_region_predictionDF$Prediction_2, " (", round(bcftools_region_predictionDF$Prob_2, digits = 5),")", sep = "")
bcftools_region_predictionDF$Region_Prediction_3 <- paste(bcftools_region_predictionDF$Prediction_3, " (", round(bcftools_region_predictionDF$Prob_3, digits = 5),")", sep = "")

bcftools_continent_predictionDF$continent_Prediction_1<- paste(bcftools_continent_predictionDF$Prediction_1, " (", round(bcftools_continent_predictionDF$Prob_1, digits = 5),")", sep = "")
bcftools_continent_predictionDF$continent_Prediction_2 <- paste(bcftools_continent_predictionDF$Prediction_2, " (", round(bcftools_continent_predictionDF$Prob_2, digits = 5),")", sep = "")
bcftools_continent_predictionDF$continent_Prediction_3 <- paste(bcftools_continent_predictionDF$Prediction_3, " (", round(bcftools_continent_predictionDF$Prob_3, digits = 5),")", sep = "")

sampleOrder <- bcftools_barcodeCSV$Sample

newVec <- c()
# for(i in sampleOrder){
#     newVar <- strsplit(i, split = "_")[[1]][1]
#     # print(newVar)
#     newVec <- c(newVec, newVar)
#     # print(newVec)
# }

#WGS samples have underscores in name, so have to adjust code
for(i in sampleOrder){
    newVar <- strsplit(i, split = "_sWGA")[[1]][1]
    print(newVar)
    newVec <- c(newVec, newVar)
    print(newVec)
}

bcftools_country_predictionDF$Sample <- newVec
bcftools_region_predictionDF$Sample <- newVec
bcftools_continent_predictionDF$Sample <- newVec



#Sample | pred1 |pred2 | pred3 | hets | missing | masked


bcftools_country_predictionDF <- bcftools_country_predictionDF[,c(15,12,13,14,8,9,10)]
bcftools_region_predictionDF <- bcftools_region_predictionDF[,c(15,12,13,14,8,9,10)]
bcftools_continent_predictionDF <- bcftools_continent_predictionDF[,c(15,12,13,14,8,9,10)]
# print(bcftools_country_predictionDF)

thirdMerge <- merge(bcftools_country_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)
fourthMerge <- merge(bcftools_region_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)
sixthMerge <- merge(bcftools_continent_predictionDF, metadataDF , by.x = "Sample", by.y = "Seq_ID_final", all.x = T)


write.table(firstMerge, paste(timestamp, time_now, "gatkCountry_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)
write.table(secondMerge, paste(timestamp, time_now, "gatkRegion_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)
write.table(thirdMerge, paste(timestamp, time_now, "bcfCountry_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)
write.table(fourthMerge, paste(timestamp,time_now, "bcfRegion_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)
write.table(fifthMerge, paste(timestamp, time_now, "gatkcontinent_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)
write.table(sixthMerge, paste(timestamp,time_now, "bcfcontinent_prediction_withMetadata.txt", sep ="_"), quote = F, sep = "\t", row.names = F)


# colnames(gatk_country_predictionDF) <- c("Country_Prediction_1", "Country_Prediction_2", "Country_Prediction_3")

# round(, digits = 4)
# analysis <- args[4]



# number <- toString(args[4])

