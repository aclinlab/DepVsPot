%% simulator

% theta will be in this order:
% 1. muXenter
% 2. sigmaXenter
% 3. muXstay
% 4. sigmaXstay

% a version to generate the negative control with zero inter-fly
% variability

function simData = simFlyDecsSigmaZero(theta, numDecsAllFlies)
% returns numFlies x 3 matrix
% dimension 2 order is: approach, stay, enter

numFlies = length(numDecsAllFlies);

Xenter = theta(1); %+randn(numFlies,1)*theta(2);
Xstay = theta(3); %+randn(numFlies,1)*theta(4);
Penter = repmat(normcdf(Xenter),numFlies,1);
Pstay = repmat(normcdf(Xstay),numFlies,1);

simData = simDecs_mex(numDecsAllFlies,Penter,Pstay);


end