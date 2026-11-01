load('JulianMelissaData.mat');
numDilutions = size(decisionsTotal,1); % a cell array with numDilutions x numVoltages cells
numVoltages = size(decisionsTotal,2);
load('integV.mat');
load('fractions.mat');

poissonNoise = false;

%% initialize parameters

% all fitted parameters are duplicated along the last dimension - to compare Tm=0
% and Tm>0

if poissonNoise
    Tms = [0 0.01 0.02 0.025 0.03];
else
    % for grid search on Tm
    % Tms = [0 0.3 0.4 0.5 0.6 0.7]; % for no sqrt() 0.4 is best
    Tms = [0, 0.1, 0.15, 0.175, 0.2, 0.225, 0.25, 0.3]; %, 0.225]; % for sqrt() %, 0.25, 0.275, 0.3, 0.4];
end
if poissonNoise
    numTmConds = length(Tms);
    Tas = 0.06^2*ones(1,numTmConds);
else
    numTmConds = length(Tms);
    Tas = 0.06*ones(1,numTmConds);
end
muSs = [0; -0.33; -0.67; -1]; %[0; 0.5; 0.5; 0.5]; %
muSs = repmat(muSs,1,numTmConds); 
sigmaSs = zeros(numVoltages,numTmConds); %sqrt( sigmaXenterIntegV_finals(:,6).^2 - sigmaXenterIntegV_finals(:,2).^2 ) / T ; %%%%%%%% HERE IS WHERE YOU ARE - SHOULD IT BE * OR ?

% Get naive params from results of ABC

S = load('naiveCombinedParams_theta_w_only.mat');
naive_mu = S.weights_all_linear' * S.thetas_all_linear;

muXenter_init = naive_mu(1);
sigmaXenter_init = naive_mu(2);
slope_muXstay_init = naive_mu(3);
slope_sigmaXstay_init = naive_mu(4);

load('Fig7G_result.mat')
% get parameters for deltaXstay habituation
mudeltaXstay_integV = mean(mudeltaXstay_finals,[2 3]);
muA_quad=fitnlm(fractions,mudeltaXstay_integV,'y~b0*x^2+b1',[-1 0.2]);
muA_qcoefs = muA_quad.Coefficients.Estimate;
mudeltaXstay_finals_mean = mean(mudeltaXstay_finals,3);
sum_deltaXstays = sum(mean(mudeltaXstay_finals,3),1); % actual sum
sum_b0_only = sum(fractions.^2*muA_qcoefs(1)); %predicted with b0 only
offsets = (sum_deltaXstays - sum_b0_only) / numDilutions; % b1



learningRate = 0.002;

% set targets for Penter, Pstay in post-training epoch
[meanPenter_target_post,meanPstay_target_post,stdPenter_target_post,stdPstay_target_post] = deal(nan(numDilutions,numVoltages));
for d=1:numDilutions
    for v=1:numVoltages
        meanPenter_target_post(d,v) = mean(decsTrainedOdorEnter{d,v}(:,6),1,"omitnan");
        meanPstay_target_post(d,v) = mean(decsTrainedOdorStay{d,v}(:,6),1,"omitnan");
        stdPenter_target_post(d,v) = std(decsTrainedOdorEnter{d,v}(:,6),0,1,"omitnan");
        stdPstay_target_post(d,v) = std(decsTrainedOdorStay{d,v}(:,6),0,1,"omitnan");
    end
end




%% Fit deltas

% we are only fitting to epoch 6
% simResultAll = cell(numDilutions,numVoltages);
[mean_simResult,std_simResult] = deal(zeros(numDilutions,numVoltages,3));
rmse = zeros(2000,numTmConds);
muSshistory = zeros(4,numTmConds,2000);
T_history = zeros(2,numTmConds,2000);
meanPenter_diff_sumAll_history = zeros(numTmConds,2000);
meanPenter_diff_byV_history = zeros(numVoltages,numTmConds,2000);


% params needs the following fields:
% -- these are to be tuned
%  Ta, additive temp
%  Tm, multiplicative temp
%  muSs, mean degree of plasticity (numVoltages x 1 array)
%  sigmaSs, stdev of plasticity (numVoltages x 1 array)
%  offsets, equal and opposite change to both Xenter and Xstay caused by shock (numVoltages x 1 array)
% -- these are fixed by the naive data
%  muXenter_pre, could be a constant, or could be numDilutions x 1
%  sigmaXenter_pre, could be a constant, or could be numDilutions x 1
%  muXstay_pre, could be parameters for a line, or could be numDilutions x1
%  sigmaXstay_pre, could be a constant, or could be numDilutions x 1
%  H, habituation factor on Xstay - a constant (applied to equation to make
%   it habituation stronger for higher conc.), or could be numDilutions x 1
% -- these are flags
%  calcNaive - default false
%  numTrials - default 1
% data needs these fields:
%  fractions, fraction of KCs activated by each dilution - numDilutions x 1
%  decisionsTotal, numDilutions x numVoltages cells array, # decs per fly

params.offsets = offsets;
params.H = muA_qcoefs(1);

params.muXenter_pre =    muXenter_init;
params.sigmaXenter_pre = sigmaXenter_init;
params.muXstay_pre_slope = slope_muXstay_init;
params.sigmaXstay_pre_slope = slope_sigmaXstay_init;

params.numTrials = 10;
params.poissonNoise = poissonNoise;

dataP.fractions = fractions;
dataP.decisionsTotal = decisionsTotal;

runningTotalNumIters = 0;
clear paramsByTmCond
tic

for TmCond = 1:numTmConds

    Ta = Tas(TmCond);
    Tm = Tms(TmCond);
    params.Ta = Ta;
    params.Tm = Tm;

    params.muSs = muSs(:,TmCond);
    params.sigmaSs = sigmaSs(:,TmCond);

    TmCond
    converged = false;
    numIters = 0;
    while ~converged 
        if numIters && ~mod(numIters,100)
            disp(numIters)
            disp(meanPenter_diff_sumAll)
            disp(meanPenter_diff_byV(2:end))
            disp(stdPenter_diff_byV(2:end))
            disp(meanPstay_diff_byV)
            % disp(params.H)
            % disp(params.offsets)
            disp([params.Ta, params.muSs'])
        end

        % run the model
        [simResult,~] = runModel(params,dataP);

        % get the average and stdev
        for d=1:numDilutions
            for v=1:numVoltages
                mean_simResult(d,v,:) = mean(simResult{d,v},[1 2],'omitnan');
                std_simResult(d,v,:) = std(simResult{d,v},0,[1 2],'omitnan');
            end
        end

        converged = true;
        % check all the target conditions. If any are not met, then set
        % converged flag back to false
        meanPenter_diff = meanPenter_target_post - mean_simResult(:,:,3); % diff>0 if you are too low % this is a matrix across all d and v
        stdPenter_diff = stdPenter_target_post - std_simResult(:,:,3);
        meanPstay_diff = meanPstay_target_post - mean_simResult(:,:,2); % diff>0 if you are too low % this is a matrix across all d and v
        % display to see the numbers hopefully doing down over GD
        % meanPenter_diff

        % update T
        meanPenter_diff_sumAll = sum(meanPenter_diff(2:end,2:end),'all','omitnan'); % diff > 0 if Penter is too low. Ignore V=0 and odor=0
        if abs(meanPenter_diff_sumAll) > 0.02
            if poissonNoise
                params.Ta = max([0, params.Ta + 0.05*learningRate*meanPenter_diff_sumAll]); % diff>0 if Penter is too low -> T should be bigger
            else
                params.Ta = max([0, params.Ta + learningRate*meanPenter_diff_sumAll]); % diff>0 if Penter is too low -> T should be bigger
            end
            converged = false;
        end

        % update muSs, sigmaSs
        meanPenter_diff_byV = sum(meanPenter_diff(2:end,:),1,'omitnan'); % sum over dilutions, ignore odor=0
        stdPenter_diff_byV = sum(stdPenter_diff(2:end,:),1,'omitnan'); % sum over dilutions, ignore odor=0

        % start v at 2 so as not to fit muS and sigmaS for 0 V, in order to
        % reduce the number of parameters that need fitting
        for v=2:numVoltages
            if abs(meanPenter_diff_byV(v)) > 0.05
                params.muSs(v) = params.muSs(v) + learningRate*meanPenter_diff_byV(v);
                converged = false;
            end

            if abs(stdPenter_diff_byV(v)) > 0.05 %0.15
                params.sigmaSs(v) = max([0 params.sigmaSs(v) + learningRate*stdPenter_diff_byV(v)]);
                converged = false;
            end

        end

        % update offsets
        meanPstay_diff_byV = sum(meanPstay_diff,1,'omitnan'); % sum over dilutions
        
        for v=1:numVoltages
            if abs(meanPstay_diff_byV(v)) > 0.05
                params.offsets(v) = params.offsets(v) + learningRate*meanPstay_diff_byV(v);
                converged = false;
            end
        end

        % update H
        sum_meanPstay_diff_all = sum(meanPstay_diff_byV,'omitnan'); % sum over voltages also

        if abs(sum_meanPstay_diff_all) > 0.05
            params.H = params.H + learningRate*sum_meanPstay_diff_all;
        end

        % % keep track of errors, let's stop when average decrease in error in
        % % the last 10 iterations is < some number
        numIters = numIters+1;
        rmse(numIters,TmCond) = mean(meanPenter_diff.^2,'all','omitnan');
        muSshistory(:,TmCond,numIters) = muSs(:,TmCond);

        T_history(:,TmCond,numIters) = [params.Ta params.Tm];
        meanPenter_diff_sumAll_history(TmCond,numIters) = meanPenter_diff_sumAll;
        meanPenter_diff_byV_history(:,TmCond,numIters) = meanPenter_diff_byV;
        % disp(params.c)
    end
    
    Tas(TmCond) = params.Ta;
    Tms(TmCond) = params.Tm;
    numIters
    runningTotalNumIters = runningTotalNumIters + numIters;
    paramsByTmCond(TmCond) = params; % store the parameters
end
runningTotalNumIters
toc

disp('done MC gradient descent')


%% run simulation with these parameters
numTrials=200;
pvalues_post = nan(numDilutions,numVoltages,3,numTmConds);
pvalues_delta = nan(numDilutions,numVoltages,numTmConds);
interSimCMD = nan(numDilutions,numVoltages,numTrials*(numTrials-1)/2,3,numTmConds);
interSimMeanCMD = nan(numDilutions,numVoltages,numTrials,3,numTmConds);
dataVsSimCMD = nan(numDilutions,numVoltages,3,numTmConds);
[simDataMC_pre,simDataMC_post] = deal(cell(numDilutions,numVoltages,numTmConds));
[mean_simDataMC, std_simDataMC, mean_simDataMC_delta, std_simDataMC_delta] = deal(nan(numDilutions,numVoltages,3,numTmConds));
[ci95_simDataMC,ci95_simDataMC_delta] = deal(nan(numDilutions,numVoltages,2,3,numTmConds));
meandeltaXenter = nan(numDilutions,numVoltages,numTmConds);
tic

for TmCond=1:numTmConds
    params = paramsByTmCond(TmCond);
    params.calcNaive = true;
    params.numTrials = numTrials;

    [simDataMC_post(:,:,TmCond),simDataMC_pre(:,:,TmCond)] = runModel(params,dataP);

    for d=1:numDilutions
        for v=1:numVoltages
            % meandeltaXenter(d,v,TmCond) = mean(deltaXenter_trials,'omitnan');
            mean_per_trial = mean(simDataMC_post{d,v,TmCond},1,'omitnan');
            mean_delta_per_trial = mean(simDataMC_pre{d,v,TmCond}-simDataMC_post{d,v,TmCond},1,'omitnan');
            mean_simDataMC(d,v,:,TmCond) = mean(mean_per_trial,2,'omitnan');
            std_simDataMC(d,v,:,TmCond) = std(mean_per_trial,0,2);
            ci95_simDataMC(d,v,:,:,TmCond) = prctile(mean_per_trial,[2.5 97.5]);

            mean_simDataMC_delta(d,v,:,TmCond) = mean(mean_delta_per_trial,2,'omitnan');
            std_simDataMC_delta(d,v,:,TmCond) = std(mean_delta_per_trial,0,2);
            ci95_simDataMC_delta(d,v,:,:,TmCond) = prctile(mean_delta_per_trial,[2.5 97.5]);
        end
    end

end
toc


%% Evaluate models: Check KS2 test with the final parameter values
for TmCond=1:numTmConds
    TmCond
    for d=1:numDilutions
        for v=1:numVoltages
            % disp([TmCond d v])
            tic
            % Use the commented out lines to swap between KS vs CvM test
            [pvalues_post(d,v,1,TmCond), interSimCMD(d,v,:,1,TmCond), dataVsSimCMD(d,v,1,TmCond), interSimMeanCMD(d,v,:,1,TmCond)] = MonteCarloCM2test(trainedOdorDecs{d,v}(:,6), simDataMC_post{d,v,TmCond}(:,:,1));
            [pvalues_post(d,v,2,TmCond), interSimCMD(d,v,:,2,TmCond), dataVsSimCMD(d,v,2,TmCond)] = MonteCarloCM2test(decsTrainedOdorStay{d,v}(:,6), simDataMC_post{d,v,TmCond}(:,:,2));
            [pvalues_post(d,v,3,TmCond), interSimCMD(d,v,:,3,TmCond), dataVsSimCMD(d,v,3,TmCond), interSimMeanCMD(d,v,:,3,TmCond)] = MonteCarloCM2test(decsTrainedOdorEnter{d,v}(:,6), simDataMC_post{d,v,TmCond}(:,:,3));
            toc

            disp(append('finished CM test for Tm condition ', num2str(TmCond), ', dilution ', num2str(d), ' voltage ', num2str(v)))
        end
    end
end

% pvalues_post
% rows: P(app),P(stay),P(enter)
% columns: Tm=0, Tm>0
%% display average divergence
dataVsSimDiv = squeeze(mean(dataVsSimCMD(2:end,2:end,:,:),[1 2]))
withinSimDiv = squeeze(mean(interSimCMD(2:end,2:end,:,:,:),[1 2 3]))
figure
subplotLabels = {'P(appr)','P(stay)','P(enter)'};
legendLabels = {'data vs. sim', 'sim vs. sim'};
for p=1:3
    subplot(1,3,p)
    title(subplotLabels{p})
    hold on
    
    plot(Tms,dataVsSimDiv(p,:),'o-','MarkerSize',3)
    plot(Tms,withinSimDiv(p,:),'o-','MarkerSize',3)
    xlabel('T_m')
    ylabel('average C-vM divergence')
    % legend(legendLabels)
    set(gca,'units','centimeters','position',[1+3*(p-1),1,2,2])
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(6,'points')
end
totalDiffDiv = sum(dataVsSimDiv-withinSimDiv,1);
[~, bestTmIndex] = min(totalDiffDiv);

%% save results
timeStamp = datestr(now,'yymmddTHHMMSS');

save(strcat('unifiedModelResults',timeStamp,'.mat'),'paramsByTmCond','dataP','bestTmIndex','pvalues_post','dataVsSimDiv','withinSimDiv',...
    'simDataMC_post','simDataMC_pre','mean_simDataMC','ci95_simDataMC') %'interSimKSD','dataVsSimKSD','interSimMeanKSD',...
