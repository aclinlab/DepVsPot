clear
load('JulianMelissaData.mat')
load('integV.mat')

options = struct;
options.N = 2000;
options.NT = options.N/2;
options.alpha = 0.95;
options.eps_target = 0.05;
options.maxIter = 2000;
options.M = 1;
options.verbose = true;
options.bisection_tol = 1e-3;
options.bisection_maxits = 40;
options.annotation = 'M=1, fixed dist bug, weighted cov, adaptive s a_target= 0.05';
options.min_stable_iters = 5;
options.accept_rate_min = 0.005; 
options.resimulate_after_resampling = true;
options.min_stable_epsilon_iters = 10;
options.rel_eps_diff_tol = 1e-4;
options.increaseM = true;
options.modelSelection = true;
options.modelCol = 5;
options.activeCols = {[1 3], [1 2 3], [1 2 3 4]}; %{1:4, 1:4, 1:4}; % use this as negative control for model selection where the 3 models are identical

load('fractions.mat')

e = 2; % naive epoch

numDilutions = 7;

naiveDecsCounts = cell(numDilutions,1);
for d=1:7
    naiveDecsCounts{d} = decisionsTotalIntegV{d}(:,e);
end

S = load('naiveCombinedSigmaZero.mat'); % neg ctrl for no sigma at all
fictParamSettings(1,1:4) = S.out.weights' * S.out.theta;
fictParamSettings(1,[2 4]) = 0;
fictParamSettings(1,5) = 1; % model 1

S = load('naiveCombinedSlopeZero.mat'); % neg ctrl for slope_sigmaXstay
fictParamSettings(2,1:4) = S.out.weights' * S.out.theta;
fictParamSettings(2,4) = 0;
fictParamSettings(2,5) = 2; % model 2

numConds = size(fictParamSettings,1) + 1;
mu = 0.2;
b = 1;

disp(datestr(now,'yymmddTHHMMSS'));

obsData = cell(numDilutions,1);
for d=1:7
    numFlies = length(naiveDecsCounts{d});
    obsData{d} = zeros(numFlies,3);
    obsData{d}(:,1) = trainedOdorDecsIntegV{d}(:,e);
    obsData{d}(:,2) = decsTrainedOdorStayIntegV{d}(:,e);
    obsData{d}(:,3) = decsTrainedOdorEnterIntegV{d}(:,e);
end

% --- Define Function Handles ---

% Kernel Sampler (The "Rejuvenation" Logic)
kernelSamplerFn = @(theta_old, cholC) utils_naive_combined_MS.kernelSamplerFn(theta_old, cholC);

% Prior Sampler
priorSamplerFn = @(n) utils_naive_combined_MS.priorSamplerNaiveFn(n, mu, b);

% Prior Log PDF
priorLogPdfFn = @(theta) utils_naive_combined_MS.priorLogPdfNaiveFn(theta, mu, b);

% % negative control simulators
% % Simulator for zero slope
% simulatorFn = @(theta) simFlyDecsNaive_slopeSigmaStayZero(theta, naiveDecsCounts, fractions);
% % Simulator for zero sigma
% simulatorFn = @(theta) simFlyDecsNaiveSigmaZero(theta, naiveDecsCounts, fractions);

% repetitions
for p = 1:10 

    for c = 1:numConds

        % % for downsampling
        % downsampleFactor = 0.5;
        % obsData = cell(numDilutions,1);
        % for d=1:7
        %     numFlies = size(decisionsTotalIntegV{d},1);
        %     downsampledNum = floor(numFlies*downsampleFactor);
        %     downsampledIndices = randsample(numFlies,downsampledNum); % sample without replacement
        %     obsData{d} = zeros(downsampledNum,3);
        %     obsData{d}(:,1) = trainedOdorDecsIntegV{d}(downsampledIndices,e);
        %     obsData{d}(:,2) = decsTrainedOdorStayIntegV{d}(downsampledIndices,e);
        %     obsData{d}(:,3) = decsTrainedOdorEnterIntegV{d}(downsampledIndices,e);
        %     naiveDecsCounts{d} = decisionsTotalIntegV{d}(downsampledIndices,e);
        % end

        % Simulator
        simulatorFn = @(theta) simFlyDecsNaive_MS(theta, naiveDecsCounts, fractions);


        if c<=size(fictParamSettings,1)
            % controls
            obsData = simFlyDecsNaive_MS(fictParamSettings(c,:), naiveDecsCounts, fractions);
        else %if c==3
            % data
            obsData = cell(numDilutions,1);
            for d=1:7
                numFlies = length(naiveDecsCounts{d});
                obsData{d} = zeros(numFlies,3);
                obsData{d}(:,1) = trainedOdorDecsIntegV{d}(:,e);
                obsData{d}(:,2) = decsTrainedOdorStayIntegV{d}(:,e);
                obsData{d}(:,3) = decsTrainedOdorEnterIntegV{d}(:,e);
            end
        end

        % Distance Function
        distanceFn = @(sim) utils_naive_combined_MS.CvMdivergenceNaive(sim, obsData);

        % order of parameters: simulator, distance, prior sampler
        result = abc_smc_adaptive_engine(simulatorFn, distanceFn, priorSamplerFn, ... %  energy_distance_on_scaled,
            priorLogPdfFn, kernelSamplerFn, options);

        out(c,p) = result;

        % save outputs on every loop in case computer crashes - delete extras after finished
        timeStamp = datestr(now,'yymmddTHHMMSS');
        save(strcat('/Users/aclin/Google Drive/Lab/Andrew_matlab/simulateMultiplex/abcoutput/abcoutput-',timeStamp,'.mat'),'out');
    end
end

beep

% can adapt the ABC diagnostics check from run_abc_sim_per_condition_MS
