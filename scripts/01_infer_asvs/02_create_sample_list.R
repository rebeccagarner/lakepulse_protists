# Create list of samples
rm(list = ls())
myfiles <- list.files(pattern = "_R1.fastq.gz")
myfiles
(samples <- gsub("_R1.fastq.gz", "", myfiles))
write.table(samples, file = "samples", quote = FALSE, row.names = FALSE, col.names = FALSE)
