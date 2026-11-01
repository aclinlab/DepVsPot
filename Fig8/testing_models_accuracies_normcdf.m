function [p_ra,p_ra_trials]= testing_models_accuracies_normcdf(Wop,classAction1,numTrials,numtrainingSamples,Y,SubDiv)
% Calculate the probability of right action (p_ra) 
% inputs:
% Wop: KC output weights (W op = Weight output)
% classAction1: which stimuli are rewarded vs punished
% numTrials, numtrainingSamples: treat the first numtrainingSamples of Y as
%  for training only, so test only on the trials after that
% Y: KC activity (binarised)
% SubDivd: 'sub' for absolute opponency, 'div' for relative opponency
% adapted by Andrew Lin from code written by Nada Abdelrahman

numm=numTrials-numtrainingSamples;

switch SubDiv
    case 'sub'
        Ta=0.0269*size(Y,1); 
        Tm=0;

    case 'div'
        Ta= 0.0215*size(Y,1); 
        Tm= 0.1710; 
    otherwise
        error('last parameter must be sub or div')
end

p_ra_trials = zeros(numm, size(Y,2));

% from FigS4E_S5.m -->  map_particles(1,[1 11])
mu_enter_0 = 1.4423 - 0.2995; % mu_enter, H0(60V)

for trials=numtrainingSamples+1:numTrials
    for odour=1:(size(Y,2))

        z1=Wop(:,1)'*Y(:,odour,trials);
        z2=Wop(:,2)'*Y(:,odour,trials);

        denom = sqrt((Ta + Tm*z1)^2 + (Ta + Tm*z2)^2);

        if (~ isempty(find(classAction1==odour,1)) )
            delta_actions=z1-z2;
        else
            delta_actions=z2-z1;
        end
        pr_action= normcdf(mu_enter_0 - delta_actions/denom);

        p_ra_trials(trials-numtrainingSamples,odour) = pr_action;
    end
end
p_ra = mean(p_ra_trials,'all');

% for debugging
if ~isfinite(p_ra)
    p_ra
    error('here')
end

end