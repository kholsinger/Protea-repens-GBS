full <- read.csv("loci_20.csv", header=TRUE, na.strings=".")
outliers <- scan("outliers-JAGS-list.txt", what="character")

include <- setdiff(colnames(full), outliers)

subset <- full[,include]

cat("Sanity checks...\n")
cat("dim(full):   ", dim(full), "\n")
cat("dim(subset): ", dim(subset), "\n")
cat("Rows equal:  ", nrow(full) == nrow(subset), "\n")
cat("Columns:     \n",
    "   full = subset+outliers",
    ncol(full) == (ncol(subset) + length(outliers)), "\n",
    "   none shared between subset & outliers",
    length(intersect(colnames(subset), outliers)) == 0, "\n")

write.csv(subset, file="no-outliers.csv", row.names=FALSE)
