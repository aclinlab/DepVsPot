function [pvalue, withinSample, CMdiv, meanWithin] = MonteCarloCM2test(data, sim)
% [pvalue, withinSample, CMdiv, meanWithin] = MonteCarloCM2test(data, sim)
% data is N x M
% sim is N x T x M
% N is the sample size
% T is how many trials were used to simulate data to compare the real data
% to
% M is the number of measurements taken from each sample. M can be 1.
% returns:
% - pvalue (M x 1) is the fraction of inter-simulation CvM divergences that are bigger
% than the CvM divergence between the true distribution and the simulations.
% That is, it's the fraction of values in "meanWithin" that are bitter than
% the mean of "vsdata".
% - withinSamples (T*(T-1)/2 x M) is the total list of T*(T-1)/2 CvM divergences within the
% simulations
% - vsdata (T x M) CvM divergences between the data vs the T simulations
% - meanWithin is the mean divergence of each simulation to all other (T-1)
% simulations

numTrials = size(sim,2);
if size(data,2) ~= size(sim,3)
    error('size(data,2) ~= size(sim,3)');
end

numMeas = size(sim,3);
pvalue = zeros(numMeas,1);

vsdata = zeros(numTrials,numMeas); % compare data vs all T trials of the simulation
withinSample = zeros(numTrials*(numTrials-1)/2, numMeas); % compare each trial to every other trial
meanWithin = zeros(numTrials,numMeas); % for each simulation trial, calculate the mean CvM test statistic between it and every other trial
for m = 1:numMeas
    within = nan(numTrials); % fill with NaNs
    for t=1:numTrials % columns in within
        [~,~,vsdata(t,m)] = cmtest2(data(:,m), squeeze(sim(:,t,m))); % Get the CvM test statistic between the data and this simulation
        for u=(t+1):numTrials % rows in within below diagonal
            [~,~,within(u,t)] = cmtest2(squeeze(sim(:,t,m)), squeeze(sim(:,u,m))); % Get the CvM test statistic between this simulation and another simulation
        end
    end
    for t=1:numTrials
        % start on column t, go from row 1 to row t-1. (if t=1, this is
        % nothing). Then start on row t, go from column t+1 to the end (if
        % you're on the last column, this is also nothing)
        meanWithin(t,m) = ( sum(within(t,1:(t-1))) + sum(within((t+1):end,t)) ) / (numTrials-1);
    end
    withinSample(:,m) = within(~isnan(within)); % only take the lower triangle which is not nans
    pvalue(m) = nnz(meanWithin(:,m) > mean(vsdata(:,m)))/size(meanWithin,1);
end
CMdiv = mean(vsdata,1); % mean across trials

end