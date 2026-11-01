function simDataAll = simFlyDecsNaive_slopeSigmaStayZero(theta, numDecsAllFlies, fractions)
% returns D x 1 cell array, each cell is a numFlies x 3 x M matrix
% dimension 2 order is: approach, stay, enter
% %%% input: 
% theta: 1 x 4
% theta will be in this order:
% 1. muXenter
% 2. sigmaXenter
% 3. slope_muXstay
% 4. slope_sigmaXstay
% assume that muXstay = -muXenter + slope_muXstay*fractions (when no odor
% is present, muXstay = -muXenter as they are opposites)
% similarly when there is no odor, sigmaXstay = sigmaXenter but sigmaXstay
% increases with increasing odor
% numDecsAllFlies: D x 1 cell array, each having all the decision counts
% fractions - fraction of KCs active at each dilution

numDilutions = length(numDecsAllFlies);
simDataAll = cell(numDilutions,1);
for d=1:numDilutions
    numFlies = length(numDecsAllFlies{d});
    simDataAll{d} = zeros(numFlies,3);
    Xenter = theta(1)+randn(numFlies,1)*theta(2);
    Xstay = (-theta(1) + theta(3)*fractions(d)) ... % muXstay
        +randn(numFlies,1)*theta(2); %sigmaXstay %(theta(2) + theta(4)*fractions(d))
    Penter = normcdf(Xenter);
    Pstay = normcdf(Xstay);
    simDataAll{d} = simDecs_mex(numDecsAllFlies{d},Penter,Pstay);
end

end