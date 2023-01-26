library(tidyverse)
library(dplyr)
library(ggplot2)

#Reading in the original dataset and the coordinates of the OCRs 
# that are present within the TADs 
OCRs_TADdomains_int <- read_delim("C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_TADdomains_int.bed", col_names = FALSE)
x_0011_df_phyloP_all <- read_csv("C:/Users/Ictinike/Documents/WrayLab/raw_data/x_0011_df_phyloP.csv")
x_0011_df_phyloP_all$chromHMM_cat_longest <- as.factor(x_0011_df_phyloP_all$chromHMM_cat_longest)

levels(x_0011_df_phyloP_all$chromHMM_cat_longest)
colnames(OCRs_TADdomains_int) <- c("seqnames", "start", "end", 
                                   "TAD_chr", "TAD_start", "TAD_end", "TAD_ID")

# Merging the two by the OCRs coordinates
OCRs_TADdomains_merged <- merge(x_0011_df_phyloP_all, OCRs_TADdomains_int, 
                                by = c("seqnames", "start", "end"))

#Need to verify that open core promoters and enhancers are within the same TAD,
# filter the dataset by these chromHMM categories
OCRs_TADdomains_prom.enhance <- OCRs_TADdomains_merged %>% 
  filter(chromHMM_cat_longest == "Active Promoter" | chromHMM_cat_longest == "Candidate Strong Enhancer")

OCRs_TADdomains_prom.enhance$TAD_ID <- as.factor(OCRs_TADdomains_prom.enhance$TAD_ID)

OCRs_TADdomains_summary <- OCRs_TADdomains_prom.enhance %>% 
  dplyr::group_by(TAD_ID) %>% 
  dplyr::summarise(Promoters = sum(chromHMM_cat_longest == "Active Promoter"),
                   Enhancers = sum(chromHMM_cat_longest == "Candidate Strong Enhancer"))

head(OCRs_TADdomains_summary)
targetTADs <- OCRs_TADdomains_summary %>% 
  filter(Promoters > 0 & Enhancers > 0)
library(data.table)
n_occur <- data.table(table(targetTADs$TAD_ID))

targetTAD_IDs <- n_occur$V1[n_occur$N > 0]

unique_TADs <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

OCRs_TADdomains_merged_targeted <- subset(OCRs_TADdomains_merged, TAD_ID %in% targetTAD_IDs)

#Remove nested TADs
OCRs_TADdomains_merged_targeted<-OCRs_TADdomains_merged_targeted[ OCRs_TADdomains_merged_targeted$TAD_ID %in% unique_TADs$TAD_ID, ]


write_csv(OCRs_TADdomains_merged_targeted, "C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_inTADs_unique.csv")
OCRs_TADdomains_summary%>%
  ggplot()+
  geom_histogram(aes(x=Enhancers), fill="red",alpha=.6)+
  geom_histogram(aes(x=Promoters), fill="#0000FF",alpha=.6)+
  theme_bw() +
  xlab(label = "Promoters and Enhancers")

OCRs_TADdomains_summary %>% 
  mutate(`Promoter-Enhacer Within TAD`=case_when(Promoters > 0 & Enhancers > 0 ~ "Yes",
                         TRUE~"No"))%>%
  group_by(`Promoter-Enhacer Within TAD`)%>%
  summarise(n=n())
  ggplot(aes(Pairs))