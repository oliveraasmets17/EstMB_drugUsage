# Scripts for analysing long-term drug usage effects

This folder contains analysis scripts used for analyzing the long-term effects of drug usage in the Estonian Microbiome Cohort [1]. The scripts used were used to produce the results for the manuscript “A hidden confounder for microbiome studies: medications used years before sample collection” [2]. Some additional functions in the scripts are used for other analysis not covered in the manuscript. 

## Expected folder structure

- RData (Help-files in rds format)
     -- Interim (Interim files for data analysis)
    -- Results (Analysis results)
- Data (Input data for the project)
- Figures (Figures for the manuscripts)

Analysis scripts are in the home folder. Occasionally, global paths are used by the scripts to indicate the location of source data. 

## Running the code

Run analysis codes
1.	DU_run_PERMANOVA.R
Script to run beta-diversity analysis across all drug drug classes. Output shown in Figure 2a and Supplementary Table 2. 
2.	DU_run_drugUsageAnalysis.R
Script to run mOTU-level univariate analysis. Different analysis are implemented: a) for identifying mOTUs associated with active drug usage by logistic regression and linear regression adjusted for BMI, gender and age. Output shown in Figure 2b and Supplementary Table 3.  b) mOTUs associated with past drug usage by logistic regression and linear regression adjusted for BMI, gender and age. Output shown in Figure 3a and Supplementary Table 7.  c) Model selection according to AIC. Output shown in Figure 3c and Supplementary Table 8. d) Sensitivity analysis for a) and b) by excluding recent antibotics users. 
3.	DU_run_machineLearningAnalysis.R 
Script to develop machine learning models for predicting drug usage based on microbiome. Output shown in Extended Data Fig. 2b and Supplementary Table 4.   
4.	DU_run_drugDrugCharacteristicsAnalysis.R
Script to analyze the effects of drug dosage on the microbiome composition using nested linear models. Output shown in Extended Data Fig. 2c.
5.	DU_run_variancePartitioningAnalysis.R
Script to characterize the effect of active and long-term drug usage on the overall variability of the fecal microbiome through multivariate analysis of the explained variance. Output shown in Figure 3d and Supplementary Table 9. 
6.	DU_run_deltaAnalysis.R
Script to assess the effect of drug initiation or discontinuation on the microbiome. Output shown in Figure 4 and Supplementary Table 11.
7.	DU_run_diseaseAnalysis.R
Script to assess the analyze associations between the mOTU abundances and prevalent diseases by linear regression adjusted for BMI, gender and age. Output shown in Extended Data Fig 4.  
Additional post-hoc analysis
8.	DU_run_deconfoundingAnalysis.R
Script to run post hoc analysis for Q1, Q2 and disease confounding to identify potential confounding factors for the univariate mOTU-drug associations as described by Forslund et al [3]. 
9.	DU_analyze_deconfounding.R
Script help summarize the  deconfounding analysis results produced by  DU_run_deconfoundingAnalysis.R.
Visualize results
10.	DU_visualize_populationDrugUsage.R
Script to produce figures 1b, 1c, 1d, Extended Data Figure 1 and Supplementary Table 1
11.	DU_visualize_activeUsage_univariate.R
Script to produce figures 2b, 2c, and Supplementary Tables 3, 7 based on the outputs of DU_run_drugUsageAnalysis.R and DU_analyze_deconfounding.R
12.	DU_visualize_machineLearningAnalysis.R 
Script to produce Extended Data Fig. 2b and Supplementary Table 4 based on the output of DU_run_machineLearningAnalysis.R
13.	DU_visualize_activeUsage_drugCharacteristics.R
Script to produce Extended Data Fig. 2c and Supplementary Table 6 based on the output of DU_run_drugDrugCharacteristicsAnalysis.R
14.	DU_visualize_activeUsage_multivariate.R
Script to produce Figure 2a and Supplementary Table 2 based on the output of DU_run_PERMANOVA.R 
15.	DU_visualize_pastUsage_univariate.R
Script to produce figures 3a, 3c, Extended Data Fig. 3, and Supplementary Table 8 based on the outputs of DU_run_drugUsageAnalysis.R and DU_analyze_deconfounding.R
16.	DU_visualize_pastUsage_multivariate.R
Script to produce Figure 3d and Supplementary Table 9 based on the output of DU_run_variancePartitioningAnalysis.R
17.	DU_visualize_deltaAnalysis_univariate.R
Script to produce Figure 4 and Supplementary Tables 10 and 11 based on the output of DU_run_deltaAnalysis.R
18.	DU_visualize_diseaseDeconfounding.R
Script to produce Extended Data Fig 4. based on the output of DU_run_diseaseAnalysis.R and DU_analyze_deconfounding.R
19.	DU_visualize_timeFromAB_nrDrugs.R
Script to produce Extended Data Fig 2a and Figure 3b

## R packages used

* tidyverse (v2.0.0)
* ppcor (v1.1)
* effsize (v0.8.1)
* compositions (v2.0.8)
* tidymodels (v1.2.0)
* vegan (v2.6.6.1)
* ggpubr (v0.6.0)
* ggthemes (v5.1.0)
* ggsci (v3.2.0)

## References
* [1] Aasmets O, Krigul KL, Lüll K, Metspalu A, Org E. Gut metagenome associations with extensive digital health data in a volunteer-based Estonian microbiome cohort. Nat Commun. 2022 Feb 15;13(1):869. doi: 10.1038/s41467-022-28464-9. PMID: 35169130; PMCID: PMC8847343.
* [2] Aasmets O, Taba N, Krigul KL, Andreson R, Org E. Long-term consequences of drug usage on the gut microbiome. medRxiv. doi: 10.1101/2024.07.17.24310548 
* [3] Forslund SK, Chakaroun R, Zimmermann-Kogadeeva M, Markó L, Aron-Wisnewsky J, Nielsen T, Moitinho-Silva L, Schmidt TSB, Falony G, Vieira-Silva S, Adriouch S, Alves RJ, Assmann K, Bastard JP, Birkner T, Caesar R, Chilloux J, Coelho LP, Fezeu L, Galleron N, Helft G, Isnard R, Ji B, Kuhn M, Le Chatelier E, Myridakis A, Olsson L, Pons N, Prifti E, Quinquis B, Roume H, Salem JE, Sokolovska N, Tremaroli V, Valles-Colomer M, Lewinter C, Søndertoft NB, Pedersen HK, Hansen TH; MetaCardis Consortium*; Gøtze JP, Køber L, Vestergaard H, Hansen T, Zucker JD, Hercberg S, Oppert JM, Letunic I, Nielsen J, Bäckhed F, Ehrlich SD, Dumas ME, Raes J, Pedersen O, Clément K, Stumvoll M, Bork P. Combinatorial, additive and dose-dependent drug-microbiome associations. Nature. 2021 Dec;600(7889):500-505. doi: 10.1038/s41586-021-04177-9. Epub 2021 Dec 8. PMID: 34880489.

Contact: oliver.aasmets@ut.ee
