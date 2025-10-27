# Hybrid-Statistical-Learning-for-COVID-19

This study presents a hybrid modeling pipeline for COVID-19 testing operations using 15,524 pediatric PCR tests. It jointly predicts test positivity, viral load (Ct values), and turnaround times (TATs), combining interpretable machine learning methods. The approach includes fairness audits, time-aware validation, and causal analysis of drive-through testing. Key findings highlight age and pandemic timing as major predictors, significant variability in Ct and TATs, and targeted operational improvements. The result is a robust, governance-integrated framework for optimizing lab workflows

All data are publicly available through the \texttt{medicaldata} R package on CRAN.  
The specific dataset analyzed is \texttt{covid\_testing}.  
Dataset DOI: \href{https://doi.org/10.32614/CRAN.package.medicaldata}{10.32614/CRAN.package.medicaldata}.  

install.packages("medicaldata")
library(medicaldata)
covid <- medicaldata::covid_testing
