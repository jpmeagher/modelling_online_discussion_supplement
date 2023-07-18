The code in this folder allows you to reproduce the results presented in
  "Modelling online discussion under circadian rhythms and superspreading 
  dynamics" by Joe Meagher and Nial Friel.

Running each of the R scripts in the /scripts folder will populate /models, 
  /evidence and /predictions folders with the required .RDS objects. The code 
  calls the [onlineMessageboardActivity](https://github.com/jpmeagher/onlineMessageboardActivity)
  package which accompanies the paper.
  
modelfit_and_evidence.R fits the branching process models to data included 
  within the package.
  
assess_predictive_performance.R computes the elpd and crps metrics for each
  model.
  
goodness_of_fit.R samples distributions of discussion sizes to assess the
 goodness-of-fit of each model to training data.
  
appendix_frequency_basis.R conducts an analysis assessing the appropriate 
  number of basis for the sinusoidal component of the intensity function
  which is included as an appendix to the manuscript.

Knitting figures_and_tables.Rmd produces all the figures and tables included 
  in the manuscript. This Rmd document is structured to follow the manuscript
  and includes all reported calculations.
