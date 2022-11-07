# BiMR
Simulation code for bidirectional Mendelian Randomization model 
This is the 1.0 version of the BiMR simulation code
Main simulation code are in the folder R/

After load all functions
Perform one simulation with bidirectional MR model:
  Sim_one(n_sim = 5, causal = "bi_infi", method = "all")

Perform one simulation with unidirectional MR model:
  Sim_one(n_sim = 5, causal = "uni", method = "all")
