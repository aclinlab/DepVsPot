# Fig 7

`simDecs_mex.mexmaca64` is for Macs. It is a compiled version of `simDecs` that runs much faster. If you don't have a Mac, compile it using:

`codegen simDecs -args {coder.typeof(1,[Inf,1]), coder.typeof(1,[Inf,1]), coder.typeof(1,[Inf,1])}`

Note that because Approximate Bayesian Computation and our behavioural models are inherently stochastic, this code will not reproduce the simulated data in the figures exactly, but the general trends will be the same.

## ABC-SMC code

There are 3 rounds of ABC-SMC. They all call `abc_smc_adaptive_engine.m`, using function handles for a simulator and the following functions which are in a `utils_[name].m` file:

* Distance
* Kernel for perturbations
* Prior sampler
* Prior log probability

1. Naive behavior at each dilution: 
	* `run_abc_sim_per_condition_MS` ('per condition' means analysing naive behavior for each dilution separately; MS stands for 'model selection'). This will fit the 7 dilutions and some controls. 
	* `simFlyDecs_MS`
	* `utils_naive_per_dilution_MS`
	* Results in: `naiveParticles.mat` (real data) and `naiveNullFits.mat` (zero variability)
	
2. Naive behavior combined across all dilutions:
	* `run_abc_all_naive_MS` Uses `naiveCombinedSigmaZero` and `naiveCombinedSlopeZero` as controls. Takes ~ 3 h to run 10 repetitions of 3 models (2 synthetic data controls and real data)
	* `simFlyDecsNaive_MS`. Variants `simFlyDecsNaive_slopeSigmaStayZero.m` and`simFlyDecsNaiveSigmaZero.m` are used to generate the negative controls above
	* `utils_naive_combined_MS`
	* Results in: `naiveCombinedParams_theta_w_only.mat`
 

3. Learned behavior:
	* `run_abc_learned` 
	* `simFlyDecsLearnedTmK`
	* `utils_learned`
	* Each ABC run takes 60-120 min on a 2021 MacBook Pro
	* Data files:
		* `learnedTmZero.mat` is the basis of the synthetic negative control - this is derived by running ABC while fixing $T_m=0$, i.e. set `theta(14)=0` in `simFlyDecsLearnedTmK.m`
		* `learnedRealData.mat` results of ABC on real data
		* `learnedNegCtrl.mat` results of ABC on the synthetic negative control
		* `learnedDownsampled.mat` results of ABC on downsampled real data (half the sample size)
	* Figures:
		* First run `FigS4E_S5.m` to read in the data and generate Fig S4E and S5
		* Keep contents in memory, then run `Fig7IJ.m`, which will produce: Fig 7I, 7J, S4F, S7, S8; and run `FigS6_param_scatter.m` for Fig S6
	* Note that the code uses `S` where the manuscript uses $\Delta w$ e.g. `muS` and `sigmaS` for the parameters governing the change in synaptic weight

		
## Gradient descent code: 

Fig 7G:

* gradient descent: `Fig7G_grad_descent.m`. Reads `naiveCombinedParams_theta_w_only.mat` as naive starting point. Takes ~ 30 min (20x rep gradient descent)
* plot `Fig7G_plot.m`. Reads `Fig7G_result.mat`

Fig 7H:

* gradient descent on the unified model: `Fig7H_unifiedModel.m`. Calls `runModel.m`. Results saved as `unifiedModelResults[timestamp].mat`
* `Fig7H_unifiedModel_loop.m`, simply runs `Fig7H_unifiedModel.m` many times
* Plot the CvM divergence with `Fig7H.m`