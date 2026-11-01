function simDataAll = simFlyDecsNaive_MS(theta, numDecsAllFlies, fractions)
% returns numDilutions x 1 cell array, each cell is a numFlies x 3 x M matrix
% dimension 2 order is: approach, stay, enter
% %%% input: 
% theta: 1 x 4
% theta will be in this order:
% 1. muXenter
% 2. sigmaXenter
% 3. slope_muXstay
% 4. slope_sigmaXstay
% 5. model number
%   1. 2 and 4 are both zero (no sigma)
%   2. 4 is zero (constant sigma for enter and stay)
%   3. full model
% assume that muXstay = -muXenter + slope_muXstay*fractions (when no odor
% is present, muXstay = -muXenter as they are opposites)
% similarly when there is no odor, sigmaXstay = sigmaXenter but sigmaXstay
% increases with increasing odor
% numDecsAllFlies: D x 1 cell array, each having all the decision counts
% fractions - fraction of KCs active at each dilution

% comment out this switch statement for negative control with identical
% models
switch theta(5)
    case 1
        theta([2 4]) = 0;
    case 2
        theta(4) = 0;
    case 3
        % nothing
    otherwise
        error('Unknown model index: %d', theta(5));
end

numDilutions = length(numDecsAllFlies);
simDataAll = cell(numDilutions,1);
for d=1:numDilutions
    numFlies = length(numDecsAllFlies{d});
    simDataAll{d} = zeros(numFlies,3);
    Xenter = theta(1)+randn(numFlies,1)*theta(2);
    Xstay = (-theta(1) + theta(3)*fractions(d)) ... % muXstay
        +randn(numFlies,1)*max(0, theta(2) + theta(4)*fractions(d)); %sigmaXstay
    Penter = normcdf(Xenter);
    Pstay = normcdf(Xstay);
    simDataAll{d} = simDecs_mex(numDecsAllFlies{d},Penter,Pstay);
end

end