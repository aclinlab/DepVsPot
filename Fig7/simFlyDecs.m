%% simulator

% theta will be in this order:
% 1. muXenter
% 2. sigmaXenter
% 3. muXstay
% 4. sigmaXstay

function simData = simFlyDecs(theta, numDecsAllFlies)
% returns numFlies x 3 matrix
% dimension 2 order is: approach, stay, enter

numFlies = length(numDecsAllFlies);

Xenter = theta(1)+randn(numFlies,1)*max(0,theta(2));
Xstay = theta(3)+randn(numFlies,1)*max(0,theta(4));
Penter = normcdf(Xenter);
Pstay = normcdf(Xstay);

simData = simDecs_mex(numDecsAllFlies,Penter,Pstay);


end