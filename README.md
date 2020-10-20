# PIMD

## “PIMD: An Integrative Approach for Drug Repositioning Using Multiple Characterization Fusion”

Abstract：The accumulation of various types of drug informatics data and computational approaches for drug repositioning can accelerate pharmaceutical research and development. However, the integration of multi-dimensional drug data for precision repositioning remains a pressing challenge. Here, we propose a systematic framework named PIMD to predict drug therapeutic properties by integrating multi-dimensional data for drug repositioning. In PIMD, drug similarity networks based on chemical, pharmacological, and clinical data are fused into an integrated drug similarity network (iDSN) composed of many clusters. Rather than simple fusion, PIMD offers a systematic way to annotate clusters. Unexpected drugs within clusters and drug pairs with a high iDSN similarity score are therefore identified to predict novel therapeutic uses. PIMD provides new insights into the universality, individuality, and complementarity of different drug properties by evaluating the contribution of each property data. To test the performance of PIMD, we use chemical, pharmacological, and clinical properties to generate an iDSN. Analyses of the contributions of each drug property indicate that this iDSN was driven by all data types and performs better than other drug similarity networks. Within the top 20 recommended drug pairs, 7drugs have been reported to be repurposed.

# Get Started

## Run Example

main.R：the main procedure of PIMD.

enrichment.R：include main enrichment analysis of PIMD.

drug properties.R：analyzes physicochemical properties of drugs.

ATC network.R：realize the  construction of ATC network.
