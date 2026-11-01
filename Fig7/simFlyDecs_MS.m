%% simulator

% theta will be in this order:
% 1. muXenter
% 2. sigmaXenter
% 3. muXstay
% 4. sigmaXstay
% 5. model number
%   1. sigmas fixed at 0
%   2. separate sigmaXenter and sigmaXstay

function simData = simFlyDecs_MS(theta, numDecsAllFlies)
% returns numFlies x 3 matrix
% dimension 2 order is: approach, stay, enter

switch theta(5)
    case 1
        theta([2 4]) = 0;
    case 2
        % nothing
    otherwise
        error('Unknown model index: %d', theta(5));
end

numFlies = length(numDecsAllFlies);

    Xenter = theta(1)+randn(numFlies,1)*max(0,theta(2));
    Xstay = theta(3)+randn(numFlies,1)*max(0,theta(4));
    Penter = normcdf(Xenter);
    Pstay = normcdf(Xstay);
    simData = simDecs_mex(numDecsAllFlies,Penter,Pstay);

end