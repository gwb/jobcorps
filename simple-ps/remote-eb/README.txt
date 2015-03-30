
# Before running the parallel simulation

1) First step is to estimate the parameters of the model, \hat{theta}, using the script:
   `Rscript run-em.R`

2) We then obtain the eb posterior distribution of tau given theta =\hat{theta}:
   `Rscript run-ref-tau.R`

3) We generate bootstrap replicates of the data:
   `Rscript bootstrap-data.R`


# The actual parallel simulation

4) From batch script, 
   `Rscript run-eb.R $integer`
