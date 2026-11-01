%% read data
clear
load('JulianMelissaData.mat');
numDilutions = size(decisionsTotal,1); % decisionsTotal is a cell array with numDilutions x numVoltages cells
numVoltages = size(decisionsTotal,2);
numEpochs = 6;

load('integV.mat')
load('fractions.mat');
S = load('naiveCombinedParams_theta_w_only.mat');
theta = S.weights_all_linear' * S.thetas_all_linear; % because there are 10 repetitions in this file

muXenter_init = theta(1);
sigmaXenter_init = theta(2);
slope_muXstay_init = theta(3);
slope_sigmaXstay_init = theta(4);

muXstay_init = NaN;
sigmaXstay_init = NaN;

%% gradient descent to find best deltaXenter etc values
% initialise variables
simDataMC_postGD4 = cell(numDilutions,numVoltages);
% simDataMC stores the simulated data (MC stands for Monte Carlo) - used to
%  evaluate the models
numRuns = 20;
[meanPenter_target,meanPstay_target,stdPenter_target,stdPstay_target] = deal(zeros(numDilutions,numVoltages));
% these are the targets for gradient descent

[mudeltaXenter_finals, mudeltaXstay_finals, sigmadeltaXenter_finals, sigmadeltaXstay_finals] = deal(zeros(numDilutions,numVoltages, numRuns));
% these store the final values of the key parameters after gradient descent
% is finished
[XenterAllSim, XstayAllSim] = deal(cell(numDilutions,numVoltages));
GDtimes = zeros(numDilutions,numVoltages, numRuns);

e = 6;

for v=1:numVoltages
    for d= 1:numDilutions
        meanPenter_target(d,v) = mean(decsTrainedOdorEnter{d,v}(:,6),1,"omitnan");
        meanPstay_target(d,v) = mean(decsTrainedOdorStay{d,v}(:,6),1,"omitnan");
        stdPenter_target(d,v) = std(decsTrainedOdorEnter{d,v}(:,6),0,1,"omitnan");
        stdPstay_target(d,v) = std(decsTrainedOdorStay{d,v}(:,6),0,1,"omitnan");
    end
end
    
parfor n=1:numRuns
    n
    for v=1:numVoltages
        for d= 1:numDilutions
            muXstay_init = -muXenter_init + slope_muXstay_init*fractions(d);
            sigmaXstay_init = sigmaXenter_init + slope_sigmaXstay_init*fractions(d);
            disp([v d])
            if ~isempty(decisionsTotal{d,v})

                numFlies = size(decisionsTotal{d,v},1);

                % initialize parameters
                
                muXenter = norminv(meanPenter_target(d,v)); % default normal cdf with mean 0, sigma 1
                muXstay = norminv(meanPstay_target(d,v));
                % we start sigmaX parameters at 0 because there's no
                % reasonable way to estimate it from the data - you can't
                % take the stdev of the norminv or logit of the data,
                % because many flies have decision ratios of 0 or 1, for
                % which the norminv/logit would be -Inf and +Inf

                % base deltaX on initial muX from the posterior of ABCSMC from
                % integV naive data

                mudeltaXenter = muXenter - theta(1);
                mudeltaXstay = muXstay - theta(3);
                sigmadeltaXenter = 0; 
                sigmadeltaXstay = 0;

                converged = false;
                learningRate = 0.01;
                tic

                % special conditions to prevent the gradient descent from getting stuck
                if v==4 && d==6 && e==2
                    convergenceLimit = 0.001;
                elseif (v==1 && d==3 && e==6)
                    convergenceLimit = 0.005;
                elseif (v==4 && d==6 && e==6) || (v==4 && d==4 && e==6)
                    convergenceLimit = 0.01;
                else
                    convergenceLimit = 0.002;
                end

                numIters = 0;
                paramRecord = []; % for debugging
                % Start the gradient descent
                while ~converged

                    Xenter_init = muXenter_init + randn(numFlies,1)*sigmaXenter_init;
                    Xstay_init = muXstay_init + randn(numFlies,1)*sigmaXstay_init;

                    % simulate numFlies random Xenter and Xstay values
                    % based on current muX and sigmaX
                    deltaXenter = mudeltaXenter + randn(numFlies,1)*sigmadeltaXenter;
                    deltaXstay = mudeltaXstay + randn(numFlies,1)*sigmadeltaXstay;
                    Xenter = Xenter_init + deltaXenter;
                    Xstay = Xstay_init + deltaXstay;

                    % convert Xenter, Xstay to P(enter), P(stay)
                    simPenter = normcdf(Xenter);
                    simPstay = normcdf(Xstay);
                    % simulate the flies
                    simResult = simDecs( decisionsTotal{d,v}(:,6), simPenter, simPstay); % columns: Papp, Pstay, Penter

                    converged = true;
                    % check all the target conditions. If any are not met, then set
                    % converged flag back to false
                    meanPenter_diff = meanPenter_target(d,v) - mean(simResult(:,3),'omitnan'); % diff>0 if you are too low
                    meanPstay_diff = meanPstay_target(d,v) - mean(simResult(:,2),'omitnan'); % diff>0 if you are too low
                    stdPenter_diff = stdPenter_target(d,v) - std(simResult(:,3),'omitnan'); % diff>0 if you are too low
                    stdPstay_diff = stdPstay_target(d,v) - std(simResult(:,2),'omitnan'); % diff>0 if you are too low

                    % update values according to gradient descent
                    if abs(meanPenter_diff) > convergenceLimit
                        mudeltaXenter = mudeltaXenter + learningRate*meanPenter_diff;% diff>0 if you are too low
                        converged = false;
                    end
                    if abs(meanPstay_diff) > convergenceLimit
                        mudeltaXstay = mudeltaXstay + learningRate*meanPstay_diff;% diff>0 if you are too low
                        converged = false;
                    end
                    if abs(stdPenter_diff) > convergenceLimit
                        newvalue = sigmadeltaXenter + learningRate*stdPenter_diff;
                        if newvalue<0
                            sigmadeltaXenter = 0; % don't let this be the thing stopping you from converging because you'll never get the sigma to be <0
                        else
                            sigmadeltaXenter = newvalue;
                            converged = false;
                        end
                    end
                    if abs(stdPstay_diff) > convergenceLimit
                        newvalue = sigmadeltaXstay + learningRate*stdPstay_diff;
                        if newvalue<0
                            sigmadeltaXstay = 0;
                        else
                            sigmadeltaXstay = newvalue;
                            converged = false;
                        end
                    end
                    % % for debugging
                    % disp([muXenter, muXstay, sigmaXenter, sigmaXstay])
                    % disp([mudeltaXenter, mudeltaXstay, sigmadeltaXstay, sigmadeltaXenter])

                    % paramRecord = [paramRecord; [muXenter, muXstay, sigmaXenter, sigmaXstay, KSdivergence]];
                end

                Xenter = Xenter_init + deltaXenter;
                Xstay = Xstay_init + deltaXstay;

                GDtimes(d,v,n) = toc;
                % disp([mudeltaXenter, mudeltaXstay, sigmadeltaXstay, sigmadeltaXenter])
                % disp('done MC gradient descent')
                mudeltaXstay_finals(d,v,n) = mudeltaXstay;
                mudeltaXenter_finals(d,v,n) = mudeltaXenter;
                sigmadeltaXstay_finals(d,v,n) = sigmadeltaXstay;
                sigmadeltaXenter_finals(d,v,n) = sigmadeltaXenter;
            end

        end
    end
end


%% save results
timestamp = datestr(now,'yymmddTHHMMSS');

save(strcat('multiplexSimResultsFromABCSMCNaive',timestamp,'.mat'),'mudeltaXenter_finals',...
    'mudeltaXstay_finals','sigmadeltaXenter_finals','sigmadeltaXstay_finals', 'offsets','muA_qcoefs');
