clear;
load('JulianMelissaData.mat')
load('integV.mat')

options = struct;
options.N = 2000;
options.NT = options.N/2;
options.alpha = 0.95;
options.eps_target = 0.01;
options.maxIter = 500;
options.M = 10;
options.verbose = true;
options.bisection_tol = 1e-3;
options.bisection_maxits = 40;
options.positiveParams = [2 4];
options.annotation = '';
options.min_stable_iters = 5;
options.accept_rate_min = 0.005; 
options.resimulate_after_resampling = true;
options.min_stable_epsilon_iters = 10;
options.rel_eps_diff_tol = 1e-4;
options.increaseM = true;
options.modelSelection = true;
options.modelCol = 5;
options.activeCols = {[1 3], [1 2 3 4]};
options

e = 2; % naive epoch

numDilutions = 7;

mu = 0.2;

% define Kernel Sampler function handle
kernelSamplerFn = @(theta_old, cholC) utils_naive_per_dilution_MS.kernelSamplerFn(theta_old, cholC);

% Prior Sampler
priorSamplerFn = @(n) utils_naive_per_dilution_MS.priorSamplerFn(n, mu);

% Prior Log PDF
priorLogPdfFn = @(theta) utils_naive_per_dilution_MS.priorLogPdfFn(theta, mu);


fictParamSettings = [1.5 nan -0.3 nan 1; % negative control for sigmas aka pos ctrl for model 1
    1.5 nan -0.3 nan 1; % do 3 reps of each control
    1.5 nan -0.3 nan 1;
    1.5 0.5 -0.3 1 2;  % positive control for sigmas aka pos ctrl for model 2
    1.5 0.5 -0.3 1 2;
    1.5 0.5 -0.3 1 2];

numConds = size(fictParamSettings,1) + numDilutions;
disp(datestr(now,'yymmddTHHMMSS'));


for c=1:numConds
    if c<=7
        % data
        d = c;
        decsCounts = decisionsTotalIntegV{d}(:,e);
        numFlies = length(decsCounts);
        obsData = zeros(numFlies,3);
        obsData(:,1) = trainedOdorDecsIntegV{d}(:,e);
        obsData(:,2) = decsTrainedOdorStayIntegV{d}(:,e);
        obsData(:,3) = decsTrainedOdorEnterIntegV{d}(:,e);
    else
        % controls
        decsCounts = decisionsTotalIntegV{7}(:,e);
        numFlies = length(decsCounts);
        obsData = simFlyDecs_MS(fictParamSettings(c-numDilutions,:), decsCounts);
    end

    % --- Define function handles depending on the conditions ---
    % Simulator
    simulatorFn = @(theta) simFlyDecs_MS(theta, decsCounts);

    % % uncomment this to make the negative control parameters
    % simulatorFn = @(theta) simFlyDecsSigmaZero(theta, decsCounts);

    % Distance Function
    distanceFn = @(sim) utils_naive_per_dilution_MS.CvMdivergence(sim, obsData);

    % order of parameters: simulator, distance, prior sampler
    result = abc_smc_adaptive_engine(simulatorFn, distanceFn, priorSamplerFn, ...
        priorLogPdfFn, kernelSamplerFn, options);

    out(c) = result;

    % save outputs on every step in case the computer crashes. Go back
    % later to delete the earlier files
    timeStamp = datestr(now,'yymmddTHHMMSS');
    save(strcat('abcoutput-',timeStamp,'.mat'),'out','options');
end

beep

%% Check ABC diagnostics 

numPlots = 10;
for d=1:7
    figure('Name',num2str(d),'Units','normalized','Position',[(d-1)/9 0 1/9 1]);

    % epsilon
    numIters = length(out(d).epsilon);
    subplot(numPlots,1,1); plot(out(d).epsilon,'-o'); title('epsilon trace'); ylabel('\epsilon'); yscale('log')
    line([0 numIters+1], [0.05 0.05])
    xlim([0 numIters+1])

    % ESS
    subplot(numPlots,1,2); plot(out(d).ESS,'-o'); title('ESS trace'); ylabel('ESS');
    xlim([0 numIters+1])
    subplot(numPlots,1,3); plot(out(d).distinctTheta_history,'-o'); title('number of unique particles'); ylabel('No. of unique particles')
    xlim([0 numIters+1])

    % change in epsilon between iterations
    subplot(numPlots,1,4); plot(out(d).rel_eps_diff_history','-o');ylabel('rel eps diff btwn iterations'), yscale('log')
    line([0 numIters+1], [0.01 0.01])
    xlim([0 numIters+1])

    % wallclock time per iteration
    subplot(numPlots,1,5);
    hold on, plot(out(d).iter_times','-o');
    xlim([0 numIters+1])

    % L (repeat samplings when too many proposals are rejected)
    subplot(numPlots,1,6);
    plot(out(d).L_history),ylabel('L')
    xlim([0 numIters+1])

    % s, called kappa in manuscript Methods. shrink kernel sampler as p_accept
    % declines
    subplot(numPlots,1,7);
    plot(out(d).s_history),ylabel('s')
    xlim([0 numIters+1])

    % acceptance rate
    subplot(numPlots,1,8);
    plot(out(d).accept),ylabel('acceptance rate')
    xlim([0 numIters+1])
    line([0 numIters+1], [0.002 0.002])
    line([0 numIters+1], [0.005 0.005])
    yscale('log')

    [N,D,iters] = size(out(d).theta_history);
    weighted_means = zeros(D,iters);
    weighted_var = zeros(D,iters);

    for i=1:iters
        for j=1:D
            weighted_means(j,i) = mean(out(d).W_history(:,i)'*out(d).theta_history(:,j,i));
            weighted_var(j,i) = var(out(d).theta_history(:,j,i),out(d).W_history(:,i));
        end
    end

    % mean theta (should plateau as the model converges)
    subplot(numPlots,1,9)
    plot(weighted_means'), ylabel('mean(\theta)')

    % var theta (should shrink as the model converges)
    subplot(numPlots,1,10)
    plot(weighted_var'),yscale('log'), ylabel('var(\theta)')

    disp(sum(out(d).iter_times))

end


