This directory includes the scripts used to test proximity assignment of the original CERES dataset against TAD data from Hi-C.

The pipeline includes scripts from 01 to 03b (I'm sorry for the poor naming).

Below is a description of each script and its purpose:
01 - (Input: OCRs_TAD_domins_int.bed, may be requested from Wray Lab JB621 DCC) Find all overlapping TADs and remove nested TADs leaving only the largest size TAD per nested TAD grouping. 
02 -  (Input: TAD_no_nested.csv, the output from 01; OCRs_TAD_domins_int.bed; may be requested from Wray Lab JB621 DCC; x_011_df_phylops.csv; phyloP and CERES information, may be requested from wray lab) Discover Element relationship with assigned Gene and create plots that may be of interest including gene size effects and decile information
02b - DEPRECATED, Use 03 instead generally (Input: TAD_no_nested.csv, the output from 01; x_011_df_phylop_TADs.csv; TAD and phyloP and CERES information, may be requested from wray lab) Examine Promoter-Enhancer gene assignment in particular  
03 - (Input: TAD_no_nested.csv, the output from 01; x_011_df_phylop_TADs.csv; TAD and phyloP and CERES information, may be requested from Wray lab) Examine Promoter-Enhancer Interaction effects on the assignment
03b - (Input: TAD_no_nested.csv, the output from 01; x_011_df_phylop_TADs.csv; TAD and phyloP and CERES information, may be requested from Wray lab; OCRs_TAD_domins_int.bed, may be requested from Wray Lab JB621 DCC) Random Forest Tests in Jameson's style on updated TAD information. I changed the input for each RF as needed so you may wish to use different variables.

Let me know if there is anything missing.
-GK
