if("iNEXT" %in% rownames(installed.packages()) == FALSE) {install.packages("iNEXT")}
library(iNEXT)
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library('tidyverse')
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
df <- read.csv(args[1], sep="\t")
print("Reading file...")
df <- df %>%
  select(clone_id,clone_seq_count)

df <- df[!duplicated(df[ , c("clone_id")]),]
print("Calculating rarefaction... This will take some minutes")
x <- iNEXT(df$clone_seq_count)

coverage <- DataInfo(df$clone_seq_count)[4]
p <-ggiNEXT(x,color.var="Order.q") +
  xlab("Number of sequences") +
  ylab("Clone diversity") +
  annotate("text", x = 3000, y =-3000, label = paste("Coverage: ",coverage), size=8) +
  coord_cartesian(ylim=c(-0,30000),clip="off")


print("Saving the plot!")
png(str_replace(args[1],".tsv","_RAREFACTION.png"),width = 1200, height = 1200,)
print(p)
dev.off()
