# Fig 6

## Explanation of raw data

`fractions.mat` has just the variable `fractions` which has the average fraction of KCs activated by each odor concentration in Fig 6C

The contents of `JulianMelissaData.mat` are:

Each fly was subjected to the following protocol:

1. No stimulus (60")
2. Vinegar vs. air (120")
3. No stimulus (60")
4. Vinegar on both sides + 12x shocks (60")
5. No stimulus (120")
6. Vinegar vs. air (120")

Each of these 6 stages is an "epoch". 

The raw behavioural data is stored in `JulianMelissaData.mat`. The variables are:

* `decisionsTotal`: the total number of decisions each fly made in each epoch
* `decsTrainedOdorEnter`: the fraction of decisions made to **enter** the side with the trained odor (out of all decisions **starting from the other side**), for each fly in each epoch
* `decstrainedOdorStay`: the fraction of decisions made to **stay on** the trained-odor side (out of all decisions **starting from that side**), for each fly in each epoch
* `trainedOdorDecs`: the fraction of decisions made in favour of the trained odor
* `speedMean`: the average walking speed (mm/s) of each fly during each epoch
* `decsStartTrained`: the total number of decisions starting on the side of the trained odor

Each is a 7x4 cell array. The rows are the different dilutions of vinegar: 0 (no odor), 1e-5, 1e-4, 1e-3, 1e-2, 1e-1 and pure vinegar. The columns are the different shock intensities: 0, 15, 30, 60 V. Each cell contains the data for that experimental condition, which is a numFlies x 6 matrix. Each row is one fly; each column is one epoch (as described above).