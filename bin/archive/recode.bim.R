# /usr/bin/Rscript
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
bim <- read.table(args[1], header = F) 
bim$V1 <- bim$V1 %>%
  dplyr::recode("I" = "1", 
                "II" = "2",
                "III" = "3",
                "IV" = "4",
                "V" = "5",
                "23" = "6",
                "MtDNA" = "9")
bim$V2 <- bim$V1
levels(bim$V2) <- c("X","I","II","III","IV","Mt","V")
bim$V2 <- paste("ce",row.names(bim),bim$V2,bim$V4,sep = "_")
write.table(bim, "gcta.recoded.bim", quote = F, col.names = F, row.names = F)