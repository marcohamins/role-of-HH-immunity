# Code for *Identifying the role of household immunity in driving individual dengue virus infection risk*

There are two major scripts, riskfit.R and modelfit.R. 
1. riskfit.R will run all analyses found in the manuscript, print important values to the console and save figures to the /Figures subfolder that were found in the manuscript.
2. modelfit.R will re-run the xgboost methods and predict infections. It will subsequently save a few data file, householdCohortData_update.rdata that can be then loaded into riskfit.R for subsequent analysis. 

Note that not all figures are reproduced due to privacy concerns. In addition, specific dates and ages have been removed from the data leading to slight variability in outcomes when implementing the training protocol. In particular the predictions found in the householdCohortData.rdata file will differ to those in the householdCohortData_update.rdata file that is created using modelfit.R. 