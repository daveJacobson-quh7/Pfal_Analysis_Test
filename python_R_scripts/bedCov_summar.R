


args <- commandArgs(trailingOnly = T)
sample <- args[1]
inFile <- read.table(args[2], sep="\t", header = F)

chroms <- c(
"AL844501.2",
"LN999944.1",
"LN999945.1",
"LN999947.1",
"AL844509.3",
"LN999946.1",
"LN999943.1",
"AL844502.2",
"AL844503.2",
"AL844504.2",
"AL844505.2",
"AL844506.3",
"AL844507.3",
"AL844508.2"
)

perChromSize <- c(640851,
947102,
1067971,
1200490,
1343557,
1418242,
1445207,
1472805,
1541735,
1687656,
2038340,
2271494,
2925236,
3291936)


total5 <- inFile[inFile[,2] >=5,3]
total20 <- inFile[inFile[,2] >=20,3]
total50 <- inFile[inFile[,2] >=50,3]
total100 <- inFile[inFile[,2] >=100,3]
total500 <- inFile[inFile[,2] >=500,3]
total1000 <- inFile[inFile[,2] >=1000,3]

outFile <- as.data.frame(matrix(ncol = 8, nrow = 14))
for(i in 1:14){
    # search <- chroms[i]
    chromLength <- perChromSize[i]
    chromFile <- inFile[grep(chroms[i], inFile[,1]),]
    totalLength <- chromFile[1,4]
    x5 <- chromFile[chromFile[,2] >=5,3]
    x20 <- chromFile[chromFile[,2] >=20,3]
    x50 <- chromFile[chromFile[,2] >=50,3]
    x100 <- chromFile[chromFile[,2] >=100,3]
    x500 <- chromFile[chromFile[,2] >=500,3]
    x1000 <- chromFile[chromFile[,2] >=1000,3]
    # print(sum(x5)/totalLength *100)
    # print(sum(x1000)/totalLength *100)
    outFile[i,1] <- chroms[i]
    outFile[i,2] <- sum(x5)/totalLength *100
    outFile[i,3] <- sum(x20)/totalLength *100
    outFile[i,4] <- sum(x50)/totalLength *100
    outFile[i,5] <- sum(x100)/totalLength *100
    outFile[i,6] <- sum(x500)/totalLength *100
    outFile[i,7] <- sum(x1000)/totalLength *100
}

# fullSize <- 23292622
# outFile[15,1] <- "FullGenome"
# outFile[15,2] <- sum(total5)/fullSize *100
# outFile[15,3] <- sum(total20)/fullSize *100
# outFile[15,4] <- sum(total50)/fullSize *100
# outFile[15,5] <- sum(total100)/fullSize *100
# outFile[15,6] <- sum(total500)/fullSize *100
# outFile[15,7] <- sum(total1000)/fullSize *100
outFile[,8] <- sample

colnames(outFile) <- c("Chromosome", "Percent_5xCoverage","Percent_20xCoverage","Percent_50xCoverage","Percent_100xCoverage","Percent_500xCoverage","Percent_1000xCoverage")

# perChromSize <- c(640851,
# 947102,
# 1067971,
# 1200490,
# 1343557,
# 1418242,
# 1445207,
# 1472805,
# 1541735,
# 1687656,
# 2038340,
# 2271494,
# 2925236,
# 3291936)





write.csv(outFile, file = paste(sample, "_3d7_perChrom_coverage.csv", sep = ""), row.names= F, quote =F)