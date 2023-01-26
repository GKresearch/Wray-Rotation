library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)

df$DHS_prop_repeat
#Read in dfs and TADs
df <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x_0011_df_phyloP_TADs.csv")
unique_TADs <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

#Remove nested TADs
df.no_overlap<-df[ df$TAD_ID %in% unique_TADs$TAD_ID, ]
df.no_overlap$CHRPOS<-paste0(df.no_overlap$seqnames,":",df.no_overlap$start,"-",df.no_overlap$end)

elements<-c("Candidate Weak Enhancer","Candidate Strong Enhancer","Active Promoter","Inactive Promoter")
elements<-df.no_overlap$chromHMM_cat_longest%>%
  unique()
element.table<-tibble(Relationship=c("Element in Multiple TADs",
                                     "Element within Gene",#"Element crosses TAD boundary",
                                     "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                     "Element-Gene Within TAD"))
##########
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
#########

test.element<-df.no_overlap%>%
  #filter(!is.na(geneStart),chromHMM_cat_longest=="Candidate Strong Enhancer")%>%
  select(c(gene_id,chromHMM_cat_longest,CHRPOS,start,end,seqnames,TAD_chr,TAD_end,TAD_start,geneStart,geneEnd,TAD_ID,distanceToTSS,geneLength))
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
library(ggpubr)
formula <- y ~ poly(x, 3, raw = TRUE)
test.element%>%
  filter(Relationship=="Element within Gene" | Relationship =="Element-Gene Within TAD")%>%
  ggplot(aes(x=geneLength/1000,y=abs(distanceToTSS)/1000,color=chromHMM_cat_longest))+
  geom_point() + 
  theme_bw() +
  geom_smooth() +
  xlab(label="Gene Length (kb)")+
  ylab(label="Distance to TSS (kb)")+
#  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")),
#                        formula = formula)+
  ggtitle(label="TAD Contains Full Element and Gene",subtitle = 
            "Gene Length to TSS Distance") +
  theme(legend.position = "none")

test.element%>%
  filter(Relationship=="Element within Gene" | Relationship =="Element-Gene Within TAD")

#t<-test.element%>%
#  filter(Relationship!="Gene off TAD")
#fwrite(test.element,"C:/Users/Ictinike/Documents/WrayLab/raw_data/x0011_Gene-Element_TAD_Domain_Check_2022-12-13.csv")
#fwrite(t,"C:/Users/Ictinike/Documents/WrayLab/raw_data/x0011_Gene-Element_TAD_Domain_Check_No-OFF-TAD_2022-12-13.csv")


#To confirm a gene with both promoter and enhancer, use the who summarise total n
#and groupby geneID and remove those geneIDs without both but only among
#Those where the TAD contains both
test.element$gene_id
#test.element<-test.element%>%
#  mutate(gene_id=geneStart+geneEnd)
gene_check <- test.element %>% 
  filter(!is.na(gene_id))%>%
  filter(Relationship!="Gene off TAD")%>%
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(Promoters = sum(chromHMM_cat_longest == "Active Promoter"),
                   Enhancers = sum(chromHMM_cat_longest == "Candidate Strong Enhancer"))

head(gene_check)
targetgenes <- gene_check %>% 
  filter(Promoters > 0 & Enhancers > 0)
library(data.table)
n_occur <- data.table(table(targetgenes$gene_id))

target_gene_IDs <- n_occur$V1[n_occur$N > 0]

test.element.iteractions<-test.element[ test.element$gene_id %in% target_gene_IDs, ]

test.element.iteractions
t<-merge(test.element.iteractions,df.no_overlap,by=c("chromHMM_cat_longest",
                                                  "start","end","seqnames",
                                                  "gene_id","distanceToTSS"))
#fwrite(t,"C:/Users/Ictinike/Documents/WrayLab/raw_data/x0011_Gene-Element_TAD_Domain_Check_Confirmed_Promoter-Enhancer_2022-12-14.csv")
test.element.iteractions
elements<-df.no_overlap$chromHMM_cat_longest%>%
  unique()
element.table<-tibble(Relationship=c("Element in Multiple TADs",
                                     "Element within Gene",#"Element crosses TAD boundary",
                                     "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                     "Element-Gene Within TAD"))
elements<-c("Candidate Weak Enhancer","Candidate Strong Enhancer","Active Promoter","Inactive Promoter")
for (el in elements){
  element.table.hold<-test.element.iteractions%>%
    filter(chromHMM_cat_longest==el)%>%
    group_by(Relationship) %>% 
    summarise(n=n()) %>% 
    mutate(percent=round((n/sum(n)*100),1)) %>%
    mutate(Element=paste0(n," (",percent,"%)"))%>%
    select(Relationship,Element)
  element.table<-merge(element.table,element.table.hold,by=c("Relationship"), all=TRUE)
}

colnames(element.table)<-c("Relationship","Cand. Enhancer (Weak)","Cand. Enhancer (Strong)","Active Promoter","Inactive Promoter")
#colnames(element.table)<-c("Relationship",elements)#("Relationship","Cand. Enhancer (Weak)","Cand. Enhancer (Strong)","Active Promoter","Inactive Promoter")
element.table$Relationship2 <- factor(element.table$Relationship, levels = c("Element within Gene","Element-Gene Within TAD",
                                                                             "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                                                             "Element in Multiple TADs"))
element.table%>%
  arrange(Relationship2)%>%
  select(-c(Relationship2))%>%
  kable()%>%
  kable_paper() #%>% 
  #pack_rows(index = c("Proximity Assignement Likely Correct" = 2, "Proximity Assignment Less Likely" = 3, "Proximity Assignment Uncertain" = 1))


test.element%>%
  ggplot(aes(x=abs(geneLength/1000),fill=Relationship))+
  geom_histogram(binwidth = 10) +
  theme_bw()+
  scale_x_continuous(limits = c(0,400))+
  #scale_y_continuous(limits = c(0,.11))+
  ylab(label="Density")+
  xlab(label="Gene Length (kb)")+
  #  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")),
  #                        formula = formula)+
  ggtitle(label="Density Plot",subtitle="Relationships against Gene Length")
