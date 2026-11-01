function [simDataMC_post, simDataMC_pre] = runModel(params, data)
% input params and data is a struct
% params needs the following fields:
% -- these are to be tuned
%  Ta, additive temp
%  Tm, multiplicative temp
%  muSs, mean degree of plasticity (numVoltages x 1 array)
%  sigmaSs, stdev of plasticity (numVoltages x 1 array)
% -- these are fixed by the naive data
%  muXenter_pre, could be a constant, or could be numDilutions x 1
%  sigmaXenter_pre, could be a constant, or could be numDilutions x 1
%  muXstay_pre, could be parameters for a line, or could be numDilutions x1
%  sigmaXstay_pre, could be a constant, or could be numDilutions x 1
% -- these are fixed from Xstay fits
%  H, habituation factor on Xstay - a constant (applied to equation to make
%   it habituation stronger for higher conc.), or could be numDilutions x 1
%  offsets, equal and opposite change to both Xenter and Xstay caused by shock (numVoltages x 1 array)
% -- these are flags
%  calcNaive - default false
%  numTrials - default 1
% data needs these fields:
%  fractions, fraction of KCs activated by each dilution - numDilutions x 1
%  decisionsTotal, numDilutions x numVoltages cells array, # decs per fly

%% initialise and check parameters
Ta = params.Ta;
Tm = params.Tm;

% flags
if isfield(params,'calcNaive') && params.calcNaive
    calcNaive = true;
else
    calcNaive = false;
end
if isfield(params,'numTrials') && params.numTrials
    numTrials = params.numTrials;
else
    numTrials = 1;
end

if isfield(params,'poissonNoise') && params.poissonNoise
    poissonNoise = true;
else
    poissonNoise = false;
end

numDilutions = size(data.decisionsTotal,1);
numVoltages = size(data.decisionsTotal,2);

fractions = data.fractions;
muXenter_pre = params.muXenter_pre;
sigmaXenter_pre = params.sigmaXenter_pre;
% the y-intercept of muXstay is -muXenter_pre because when there's no odor,
% P(stay) = 1 - P(enter), thus Xstay = -Xenter
muXstays_pre = fractions*params.muXstay_pre_slope - muXenter_pre; 
sigmaXstays_pre = sigmaXenter_pre + fractions*params.sigmaXstay_pre_slope;

Hs = params.H*fractions.^2;

muSs = params.muSs;
sigmaSs = params.sigmaSs;
offsets = params.offsets;
decisionsTotal = data.decisionsTotal;
simDataMC_post = cell(numDilutions, numVoltages);
simDataMC_pre = cell(numDilutions, numVoltages);


%% run the simulation
for d=1:numDilutions
    habXstay = Hs(d);
    fraction = fractions(d);
    % muXenter_pre = muXenters_pre(d);
    % sigmaXenter_pre = sigmaXenters_pre(d);
    muXstay_pre = muXstays_pre(d);
    sigmaXstay_pre = sigmaXstays_pre(d);
    for v=1:numVoltages
        numFlies = size(decisionsTotal{d,v},1);
        offset = offsets(v);
        % this is to store the simulation using the final parameters from the
        % gradient descent
        simDataMC_post{d,v} = zeros(numFlies,numTrials,3); % MC = Monte Carlo
        % last dimension is: 1. trainedOdorDecs. 2. decsStay. 3. decsEnter

        % calculate numTrials reps with the "final" parameter values
        for t=1:numTrials
            S = muSs(v) + randn(numFlies,1)*sigmaSs(v);
            % cap abs(S) <= 1
            S(S>1) = 1;
            S(S<-1) = -1;

            if poissonNoise
                deltaXenter = S*fraction./sqrt( 2*Ta + Tm*(2+S)*fraction ) - offset;
            else
                % alternative denominator if not explicitly modeling noise: Ta+Tm*(2+S)*fraction
                deltaXenter = S*fraction./sqrt( (Ta + Tm*fraction)^2 + (Ta + Tm*(1+S)*fraction).^2 ) - offset;
            end
            % randomly sample the naive Xenter
            Xenter_pre = muXenter_pre + randn(numFlies,1)*sigmaXenter_pre;
            % simulate post-learning Xenter
            Xenter = Xenter_pre + deltaXenter;

            Xstay_pre = muXstay_pre + randn(numFlies,1)*sigmaXstay_pre; 
            Xstay = Xstay_pre + habXstay + offset; 

            simPenter_pre = normcdf(Xenter_pre);
            simPstay_pre = normcdf(Xstay_pre);
            simPenter = normcdf(Xenter);
            simPstay = normcdf(Xstay);

            if calcNaive
                simDataMC_pre{d,v}(:,t,:) = simDecs_mex( decisionsTotal{d,v}(:,2), simPenter_pre, simPstay_pre);
            end
            simDataMC_post{d,v}(:,t,:) = simDecs_mex( decisionsTotal{d,v}(:,6), simPenter, simPstay);
        end
    end
end

end