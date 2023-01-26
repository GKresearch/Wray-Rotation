
#read in data
df <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/x_0011_df_phyloP.csv")

#make parent dataframes and summarize
ForestData <- df %>% 
  dplyr::select(wgCERES_score_nosig, signalValue, pValue, DHS_prop_repeat, 
                DHS_prop_GC, DHS_length, n_SNV_Zhou_per_bp, dependency, 
                distanceToTSS, zeta.human, zeta.chimp, PP_con, PP_acc, 
                PhastCons, ploidyZhou, probIntolerantLoF, numTKOHits_Hart, 
                H3K27ac_CPM_per_1kbp, chromHMM_cat_longest, DNase_CPM_per_1kbp, 
                annotation, dhs_0_1_wg, PhyloP_placental_score,PhyloP_mammals_score,
                PhyloP_primates_score) %>% 
  na.omit()

summary(ForestData)

#ID sites in f/u screen
Fu_data <- df %>% 
  dplyr::select(wgCERES_score_nosig, signalValue, pValue, DHS_prop_repeat, 
                DHS_prop_GC, DHS_length, n_SNV_Zhou_per_bp, dependency, 
                distanceToTSS, zeta.human, zeta.chimp, PP_con, PP_acc, 
                PhastCons, ploidyZhou, probIntolerantLoF, numTKOHits_Hart, 
                H3K27ac_CPM_per_1kbp, chromHMM_cat_longest, DNase_CPM_per_1kbp, 
                annotation, dhs_0_1_wg, PhyloP_placental_score, 
                dhs_0_1_distal) %>% 
  na.omit()
# Fu_names <- Fu_data$name %>% as.vector()
# #filter on sites that successfully validated
# Forest_fu_val <- Fu_data %>% filter(dhs_0_1_wg == dhs_0_1_distal)
# Fu_val_names <- Forest_fu_val$name %>% as.vector()

#make training response vectors 
ForestOutput_CERES <- ForestData$wgCERES_score_nosig %>% as.vector() 
ForestOutput_Signal <- ForestData$signalValue %>% as.vector()
ForestOutput_pVal <- ForestData$pValue %>% as.vector()
ForestOutput_dhs_wg_sig <- ForestData$dhs_0_1_wg %>% as.vector()

#make training input dataframe 
ForestInput <- ForestData %>% 
  dplyr::select(-c(wgCERES_score_nosig, signalValue, pValue, dhs_0_1_wg))

#testing outputs   -getting NA errors in 
FuOutput_CERES <- Fu_data$wgCERES_score_nosig %>% as.vector() 
FuOutput_Signal <- Fu_data$signalValue %>% as.vector()
FuOutput_pVal <- Fu_data$pValue %>% as.vector()
FuOutput_dhs_wg_sig <- Fu_data$dhs_0_1_wg %>% as.vector()

#testing input
FuInput <- Fu_data %>% 
  dplyr::select(-c(wgCERES_score_nosig, signalValue, pValue, dhs_0_1_wg, dhs_0_1_distal))



tic()
CERES_Forest <- randomForest(x = ForestInput, y = ForestOutput_CERES,  ntree = 500, keep.forest = T, importance = TRUE)
toc()
print(CERES_Forest)





save(CERES_Forest, file = "C:/Users/Ictinike/Documents/WrayLab/CERES_Forest.Rda")
tic()
CERES_Forest_tested <- randomForest(x = ForestInput, y = ForestOutput_CERES,  ntree = 500, xtest = FuInput, ytest = FuOutput_CERES, keep.forest = T, importance = TRUE)
toc() #expect ~3 hours
save(CERES_Forest_tested, file = "C:/Users/Ictinike/Documents/WrayLab/CERES_Forest_tested.Rda")
tic()
dhsWGsig_Forest <- randomForest(x = ForestInput, y = ForestOutput_dhs_wg_sig, ntree = 500, keep.forest = T, importance = TRUE)
toc() #expect ~3 hours
save(dhsWGsig_Forest, file = "C:/Users/Ictinike/Documents/WrayLab/dhsWGsig_Forest.Rda")
dhsWGsig_Forest_tested <- randomForest(x = ForestInput, y = ForestOutput_dhs_wg_sig, ntree = 500, xtest = FuInput, ytest = FuOutput_dhs_wg_sig, keep.forest = T, importance = TRUE)
save(dhsWGsig_Forest_tested, file = "C:/Users/Ictinike/Documents/WrayLab/dhsWGsig_Forest_tested.Rda")


pVal_Forest <- randomForest(x = ForestInput,  y = ForestOutput_pVal, ntree = 500, keep.forest = T, importance = TRUE)
save(pVal_Forest, file = "C:/Users/Ictinike/Documents/WrayLab/pVal_Forest.Rda")
pVal_Forest_tested <- randomForest(x = ForestInput, y = ForestOutput_pVal, ntree = 500, xtest = FuInput, ytest = FuOutput_pVal, keep.forest = TRUE, importance = TRUE)
save(pVal_Forest_tested, file = "C:/Users/Ictinike/Documents/WrayLab/Val_Forest_tested.Rda")
Signal_Forest <- randomForest(x = ForestInput, y = ForestOutput_Signal, ntree = 500, keep.forest = TRUE, importance = TRUE)
save(Signal_Forest, file = "C:/Users/Ictinike/Documents/WrayLab/Signal_Forest.Rda")
Signal_Forest_tested <- randomForest(x = ForestInput, y = ForestOutput_Signal, ntree = 500, xtest = FuInput, ytest = FuOutput_Signal, keep.forest = TRUE, importance = TRUE)
save(Signal_Forest_tested, file = "C:/Users/Ictinike/Documents/WrayLab/Signal_Forest_tested.Rda")

