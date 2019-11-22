# FLOMICS_Anjali

Currently code uploaded by KI for processing variant data provided by BC genomics team. This involves mainly the following:
1. Extract for each patient, a VCF file with merged variants (SNPs) from Platypus and loFreq
2. Extract for each patient, a VCF file with merged variants (INDELs) from Platypus and loFreq
3. Normalize each patient's VCF file 
4. Annotated each patient's VCF file with ANNOVAR (cosmic and population variants in addition to gene localization)
5. Combine all VCF files from patients into one matrix for further analysis 
6. Basic summary plots of all variants (after some filtering such as removing population variants)
