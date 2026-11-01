function simDataAll = simFlyDecsLearnedTmK(theta, numDecsAllFlies, fractions)
% simFlyDecsLearned
% Simulates fly decisions under electric shock conditioning using the
% multiplicative noise and Dilution-dependent habituation.
%
% Inputs:
%   theta: 1xD vector of parameters (Single particle)
%          1-4:   [muXenter_pre, sigmaXenter_pre, slope_muXstay, slope_sigmaXstay]
%          5-7:   muSs (15V, 30V, 60V)
%          8-11:  offsets (0V, 15V, 30V, 60V)
%          12:    H
%          13:    Ta
%          14:    Tm
%          15-17: sigmaSs (15V, 30V, 60V0
%          18:    k - fraction of variance in Xenter and Xstay that is not consistent between trials ("non-personality")

%   numDecsAllFlies: Data structure used to get N flies per condition and
%                    num. decisions per fly
%   fractions: fractions of active KCs
%
% Output:
%   simDataAll: {dilution,voltage} cell array of (NumFlies x 3 x 2) matrices
%               Columns: [Approach, Stay, Enter]
%               3rd dimension: pre, post

% --- unpack parameters ---
% naive
muXenter_pre     = theta(1);
sigmaXenter_pre  = theta(2);
slope_muXstay    = theta(3);
slope_sigmaXstay = theta(4);

% learned
muSs    = [0 theta(5:7)]; % muS for 0 V is forced to be 0
offsets = theta(8:11);
H = theta(12);
Ta        = theta(13); %0.09; %
Tm    = theta(14);
sigmaSs  = [0 theta(15:17)]; %[0 0 0 0 ]; %

k = theta(18); % fraction of variance in Xenter and Xstay that is not consistent between trials ("non-personality")

[numDilutions, numVoltages] = size(numDecsAllFlies);
simDataAll = cell(numDilutions, numVoltages);

for d = 1:numDilutions
    frac = fractions(d);

    % --- calculate naive parameters ---
    % mu_stay = -mu_enter + slope * n
    muXstay_pre = -muXenter_pre + slope_muXstay * frac;

    % sigma_stay = sigma_enter + slope * n (Updated per instruction)
    sigmaXstay_pre = sigmaXenter_pre + slope_sigmaXstay * frac;

    % --- calculate non-associative habituation ---
    % Quadratic function of # active KCs
    habXstay = H * frac^2; 

    for v = 1:numVoltages
        % Get numFlies from observed data to match sample size
        if isempty(numDecsAllFlies{d,v})
            continue;
        end
        numFlies = size(numDecsAllFlies{d,v}, 1);

        S = muSs(v) + randn(numFlies,1)*sigmaSs(v);
        % restrict S to be between -1 and 1
        S(S>1) = 1;
        S(S<-1) = -1;

        offset = offsets(v);

        % --- calculate associative learning (delta Xenter) ---

        if frac == 0
            deltaXenter = -offset; % No odor association possible, only offset
        else
            Avoid = max([0,Ta + Tm*frac]);
            Approach = Ta + Tm.*(1+S).*frac;
            Approach(Approach<0) = 0;
            term_denom = sqrt(Avoid^2 + Approach.^2);
            term_num   = S * frac; 
            deltaXenter = (term_num ./ term_denom) - offset;
        end

        % --- simulation ---

        % sample naive variables

        Xenter_pers = muXenter_pre + randn(numFlies, 1) * sigmaXenter_pre * sqrt(1-k);
        Xstay_pers  = muXstay_pre  + randn(numFlies, 1) * sigmaXstay_pre * sqrt(1-k);

        Xenter_pre = Xenter_pers + randn(numFlies, 1) * sigmaXenter_pre * sqrt(k);
        Xstay_pre  = Xstay_pers  + randn(numFlies, 1) * sigmaXstay_pre * sqrt(k);

        % apply learning
        Xenter_post = Xenter_pers + deltaXenter      + randn(numFlies, 1) * sigmaXenter_pre * sqrt(k);
        Xstay_post  = Xstay_pers + habXstay + offset + randn(numFlies, 1) * sigmaXstay_pre * sqrt(k);

        % convert X to P (probabilities) via normcdf
        Pstay_post  = normcdf(Xstay_post);
        Pstay_pre  = normcdf(Xstay_pre);
        Penter_post = normcdf(Xenter_post);
        Penter_pre = normcdf(Xenter_pre);

        % run the simulation
        % pre
        simDataAll{d,v}(:,:,1) = simDecs_mex(numDecsAllFlies{d,v}(:,1),Penter_pre,Pstay_pre);
        % post
        simDataAll{d,v}(:,:,2) = simDecs_mex(numDecsAllFlies{d,v}(:,2),Penter_post,Pstay_post);
    end
end
end