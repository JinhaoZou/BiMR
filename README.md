# BiMR
Simulation code for bidirectional Mendelian Randomization model         
This is the 1.0 version of the BiMR simulation code          
Simulation code are in the folder R/            
You can install in by using code:         
`library(devtools)`     
`install_github("JinhaoZou/BiMR")`  

## Overview
Other required packages can be installed through      
`install.packages("MendelianRandomization")`     
`install.packages("ivmodel")`      
Before run the code, packages should be load             
`library(MendelianRandomization)`       
`library(ivmodel)`        

There are major three functions: 
- Data.R, a function to generate simulated data
    -  Example: `dat = Data(causal = "uni")`, generating data with unidirectional causal effects between two phenotypes 
    -  Example: `dat = Data(causal = "bi_infi")`, generating data with bidirectional causal effects between two phenotypes
- Est_all.R, a funtion to estimate causal effects from generated data sets
    -  Example: `est.one = Est_all(data = dat$data, method = "all")` generating data estimation with all estimation methods
    -  method can set as "Ratio", "BiRatio", "Liml", "BiLiml" or "all" to estimate causal effects using one or all methods
- Sim_all.R, a function to generate 1000 simulation replicates with estimated causal effects for selected scenario setting 
    -  Example: `Sim_all(causal = "uni", method = "all")`

## Examples for using current method to estimate the causal effect of your interested observational individual data
The sample data is saved in folder data/



