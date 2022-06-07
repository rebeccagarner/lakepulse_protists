# Modify raw fastq file names
rm(list = ls())
(myfiles <- list.files(pattern = "fastq.gz"))
gsub(".*LP-watercolumn-18S_", "LP-watercolumn-18S_", myfiles)
file.rename(from = myfiles, to = gsub(".*LP-watercolumn-18S_", "LP-watercolumn-18S_", myfiles))
list.files()
