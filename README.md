# GibbsSampler
A project where I code my own Gibbs sampler to deepen my understanding of Bayesian approaches

## Data
Data available [here](https://www.kaggle.com/datasets/rajyellow46/wine-quality).

## Files
The files in this folder are:

* JudithNeve_report - the report (source Rmd & knitted PDF)
* JudithNeve_code - the R code necessary for the report

## Code

Methods:
* prior parameters:
	- run lines 10-15 to get the data
	- run line 17 for the prior model
	- run lines 134-137 to extract the prior estimates (Table 1)

Results:
* Data:
	- run lines 27-36
* Models 1 and 2:
	- run lines 41-120 to get the MCMC sampler function
	- run lines 123-126 to set the seed and number of iterations
	- run lines 130-141 for Model 1
	- run lines 145-156 for Model 2
* Convergence (not discussed extensively in the report, but mentioned):
	- run line 4 to load the scales package
	- model 1: run lines 171-233
	- model 2: run lines 234-294
* Model 3:
	- run lines 41-120 to get the MCMC sampler function
	- run lines 123-126 to set the seed and number of iterations
	- run lines 160-167 to sample the parameters
* Model 3 - convergence:
	- run line 4 to load the scales package
	- run lines 299-308 to obtain the traceplots (Figure 1)
	- run lines 312-328 to obtain the autocorrelation plots (Figure 2)
	- run lines 332-348 to obtain the Gelman-Rubin statistics and the maximum deviation from 1
	- run lines 352-363 to obtain the MC errors and the values they are compared to
	- run line 427 and look at Accepted.MH to get the acceptance rate of the MH algorithm
* Posterior predictive check:
	- run line 123 to set the seed
	- run lines 372-392 to obtain the simulated datasets
	- run lines 397-414 to obtain the test statistics of the simulated datasets and the real data
	- run lines 419-424 to obtain the ppp-value
* Parameter estimates:
	- rune line 361 to combine the two chains into one matrix
	- run lines 429-430 to obtain parameter estimates and credible intervals
* DIC:
	- run lines 437-463 to make the function that calculates the DIC
	- run line 466 for the DIC of Model 3
	- run line 467 for the DIC of Model 1
	- run line 468 for the DIC of Model 2
* Bayes Factor:
	- run line 3 to load the bain package
	- run line 123 to set the seed
	- run lines 472-476 to get the Bayes Factor
