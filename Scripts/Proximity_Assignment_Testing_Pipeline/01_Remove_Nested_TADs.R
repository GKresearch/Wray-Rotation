library(tidyverse)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

#Read in OCRs_TADdomains then remove all but TAD information
df <- read_delim("C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_TADdomains_int.bed", col_names = FALSE)
colnames(df) <- c("seqnames", "start", "end", 
                                   "TAD_chr", "TAD_start", "TAD_end", "TAD_ID")
df<-df%>%
  select(TAD_chr,TAD_start,TAD_end,TAD_ID)%>%
  unique()

# Set up empty GRanges dataframe
gr2 = with(df[2,], GRanges(TAD_chr, IRanges(start = TAD_start, end = TAD_end, names = TAD_ID)))
gr1 = with(df[1,], GRanges(TAD_chr, IRanges(start = TAD_start, end = TAD_end, names = TAD_ID)))
overlap.check = findOverlaps(query = gr1, subject = gr2, type = "within")
overlap.df = data.frame(gr1[queryHits(overlap.check),], gr1[subjectHits(overlap.check),])

#Confirm the overlap.df is empty
if (dim(overlap.df)[1] == 0) {
  print("Empty")
}

#This for loop should keep only the largest TAD in each nested TAD group
for(i in 1:nrow(df)){

  # Make GRanges object for each TAD against every other TAD
  gr1 = with(df[i,], GRanges(TAD_chr, IRanges(start = TAD_start, end = TAD_end, names = TAD_ID)))
  gr2 = with(df[-c(i),], GRanges(TAD_chr, IRanges(start = TAD_start, end = TAD_end, names = TAD_ID)))
  range_check = findOverlaps(query = gr1, subject = gr2, type = "within")
  range_df = data.frame(df[i,][queryHits(range_check),], df[-c(i),][subjectHits(range_check),])
  
  #If no overlap is present skip, otherwise add to the overlap.df
  if (dim(range_df)[1] == 0) {
    next
  }
  else{
    overlap.df<-rbind(overlap.df,range_df)
  }
}

#Remove any nested TADs
TAD_no_overlap<-df[ !df$TAD_ID %in% overlap.df$TAD_ID, ]

unique_TADs <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

#Remove nested TADs
df.no_overlap<-df[ df$TAD_ID %in% unique_TADs$TAD_ID, ]

#Write to file
fwrite(TAD_no_overlap,"C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")
  