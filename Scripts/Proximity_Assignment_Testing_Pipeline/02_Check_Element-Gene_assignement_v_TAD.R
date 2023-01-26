library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)

#Read in dfs and TADs
df <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x_0011_df_phyloP_TADs.csv")
unique_TADs <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

#Remove nested TADs
df.no_overlap<-df[ df$TAD_ID %in% unique_TADs$TAD_ID, ]
df.no_overlap$CHRPOS<-paste0(df.no_overlap$seqnames,":",df.no_overlap$start,"-",df.no_overlap$end)

elements<-c("Candidate Weak Enhancer","Candidate Strong Enhancer","Active Promoter","Inactive Promoter")
element.table<-tibble(Relationship=c("Element in Multiple TADs",
                                     "Element within Gene",#"Element crosses TAD boundary",
                                     "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                     "Element-Gene Within TAD"))

for (el in elements){
  test.element<-df.no_overlap%>%
    filter(!is.na(geneStart),chromHMM_cat_longest==el)%>%
    select(c(chromHMM_cat_longest,CHRPOS,start,end,seqnames,TAD_chr,TAD_end,TAD_start,geneStart,geneEnd,TAD_ID,distanceToTSS,geneLength))
  n_occur <- data.frame(table(test.element$CHRPOS))
  test.element<-test.element%>%
    mutate(element_check=case_when(CHRPOS %in% n_occur$Var1[n_occur$Freq > 1]~"Element in Multiple TADs",
                                   start<TAD_start|end>TAD_end~"Element crosses TAD boundary",
                                   TRUE~"No"))%>%
    distinct(CHRPOS,.keep_all=TRUE)%>%
    mutate(gene_check=case_when(geneStart<start&start<geneEnd~"In",
                                geneStart<end&end<geneEnd~"In",
                                TRUE~"Out"),
           Relationship=case_when(element_check=="Element in Multiple TADs"~"Element in Multiple TADs",
                                  gene_check=="In"~"Element within Gene",
                                  #element_check=="Element crosses TAD boundary"~"Element crosses TAD boundary",
                                  geneEnd<TAD_start~"Gene off TAD",
                                  geneStart>TAD_end~"Gene off TAD",
                                  geneStart<TAD_start~"Gene Start off TAD",
                                  geneEnd>TAD_end~"Gene End off TAD",
                                  TRUE~"Element-Gene Within TAD"))
  element.table.hold<-test.element%>%
    group_by(Relationship) %>% 
    summarise(n=n()) %>% 
    mutate(percent=round((n/sum(n)*100),1)) %>%
    mutate(Element=paste0(n," (",percent,"%)"))%>%
    select(Relationship,Element)
  element.table<-merge(element.table,element.table.hold,by=c("Relationship"), all=TRUE)
}

colnames(element.table)<-c("Relationship","Cand. Enhancer (Weak)","Cand. Enhancer (Strong)","Active Promoter","Inactive Promoter")
element.table$Relationship2 <- factor(element.table$Relationship, levels = c("Element within Gene","Element-Gene Within TAD",
                                                                             "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                                                             "Element in Multiple TADs"))
element.table%>%
  arrange(Relationship2)%>%
  select(-c(Relationship2))%>%
  kable()%>%
  kable_paper() %>% 
  pack_rows(index = c("Proximity Assignement Likely Correct" = 2, "Proximity Assignment Less Likely" = 3, "Proximity Assignment Uncertain" = 1))
  

test.element<-df.no_overlap%>%
  filter(!is.na(geneStart),chromHMM_cat_longest=="Active Promoter")%>%
  select(c(chromHMM_cat_longest,CHRPOS,start,end,seqnames,TAD_chr,TAD_end,TAD_start,geneStart,geneEnd,TAD_ID,distanceToTSS,geneLength))
n_occur <- data.frame(table(test.element$CHRPOS))
test.element<-test.element%>%
  mutate(element_check=case_when(CHRPOS %in% n_occur$Var1[n_occur$Freq > 1]~"Element in Multiple TADs",
                                 start<TAD_start|end>TAD_end~"Element crosses TAD boundary",
                                 TRUE~"No"))%>%
  distinct(CHRPOS,.keep_all=TRUE)%>%
  mutate(gene_check=case_when(geneStart<start&start<geneEnd~"In",
                              geneStart<end&end<geneEnd~"In",
                              TRUE~"Out"),
         Relationship=case_when(element_check=="Element in Multiple TADs"~"Element in Multiple TADs",
                                gene_check=="In"~"Element within Gene",
                                #element_check=="Element crosses TAD boundary"~"Element crosses TAD boundary",
                                geneEnd<TAD_start~"Gene off TAD",
                                geneStart>TAD_end~"Gene off TAD",
                                geneStart<TAD_start~"Gene Start off TAD",
                                geneEnd>TAD_end~"Gene End off TAD",
                                TRUE~"Element-Gene Within TAD"))
test.element%>%
  ggplot()+
  geom_point(aes(x=geneLength/1000,y=abs(distanceToTSS)/1000,color=Relationship)) + 
  theme_bw() +
  geom_smooth(aes(x=geneLength/1000,y=abs(distanceToTSS)/1000,color=Relationship)) +
  xlab(label="Gene Length (kb)")+
  ylab(label="Distance to TSS (kb)")+
  ggtitle(label="Candidate Strong Enhancer",subtitle = 
  "Gene Length to TSS Distance")


max(test.element$geneLength)
test.element.dec<-test.element%>%
  mutate(decile_length = ntile(geneLength,10))%>%
  mutate(decile_length_TSS = ntile(abs(distanceToTSS),10))

test.element.dec%>%
  filter(Relationship!="Gene off TAD")%>%
  group_by(decile_length_TSS) %>% 
  rename(`TSS - Decile`=decile_length_TSS)%>%
  summarise(N=n()) %>% 
  mutate(Percent=round((N/sum(N)*100),1))%>%
  kable()%>%
  kable_minimal()
  
test.element.dec%>%
  filter(Relationship=="Gene off TAD")%>%
  group_by(decile_length) %>% 
  rename(`Gene Length - Decile`=decile_length)%>%
  summarise(N=n()) %>% 
  mutate(Percent=round((N/sum(N)*100),1))%>%
  kable()%>%
  kable_minimal()

t2<-test.element.dec%>%
  filter(decile_length_TSS==10)
max(abs(t2$distanceToTSS))
min(abs(t2$distanceToTSS))

19711/1000
1512398/1000

test.element.dec%>%
  filter(Relationship=="Gene off TAD")%>%
  ggplot(aes(x=abs(distanceToTSS/1000)))+
  geom_histogram(fill="blue",alpha=.7) +
  xlab(label="Distance to TSS (kb)") +
  ylab(label="Count") +
  theme_bw() +
  ggtitle(label = "Gene off TAD Distance Histogram")
