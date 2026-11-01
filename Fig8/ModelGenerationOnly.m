%% this code investigates the effect of varrying network parameters
% one at a time, PN-KC weights, KCs spiking thresholds and number of inputs
% on the memory performance of model fruit-flies Drosophila
% this simulation uses realistic odor responses recordings from
% (Hallem-Carlson 2006) and transformed from ORN to PN responses as by
% Olsen et.al; detailed in the Methods section in our paper.
% This code was adapted by Andrew Lin from code written by Philippe Fischer
% for Bergmann G et al 2026 J Physiol, which was adapted from code written
% by Nada Abdelrahman for Abdelrahman et al 2021 PNAS

clear all;
% param. initializations
numModelInstances = 20; % trials, each corresponding to different mushroom body instance

numTrials = 30; %number of noisy trials per each odors
noise = 1;

%loading Hallem-Carlson data
[PNactivity,~,~] = propagateORNs2PNs;
PNactivity = PNactivity';
K_odors = size(PNactivity,2);
nPNs = size(PNactivity,1);

allTestedOdors = cell(numModelInstances,1);
MBmodels_all = cell(4, numModelInstances);

% the odors we used in Fig2, 3 and 4, stored in this .mat file.
% load('Data_submitted_fly_wNoise11.mat','PNtrials');

%% loops over different number of input odors, N=50,100,150,200.
% testing models at different number of odors, Fig2C1

% random fly instantiations
for randomTrials=1:numModelInstances
    
    % Change from Philippe's code in Bergmann G et al 2026: draw the variable PN-KC connections
    % first so they are the same for all the variable-N models in this
    % iteration

    % ->varargin for mean, std and lower/upper bounds of Nclaws and theta
    MBmodel_varN_varW_varTh = MBmodelBuilder(2000,24, true, true, true); 
    % this line is commented out because the homogeneous model isn't in the
    % manuscript but it gives basically the same results as the other
    % models
    % MBmodel_fixN_fixW_fixTh = MBmodelBuilder(2000,24, false,false,false);
    MBmodel_varN_varW_varTh_actIndepHomeo = MBmodel_varN_varW_varTh;
    MBmodel_varN_varW_varTh_actIndepHomeo.PNtoKC = drawWeights_PWgivenNtheta(MBmodel_varN_varW_varTh_actIndepHomeo);  % redraw weights according to P(W|N,theta)
    
    PNtrials = zeros(nPNs, K_odors, numTrials);
    PNtrials(:,:,1) = PNactivity; %first trial is "simply" the first instance of the fictitious odor response
    % for all trials, add some noise to the generated fictitious odor response, and calculate a SNR ratio
    for t = 1:numTrials-1
        PNtrials(:,:,t+1) = PNactivity + ... % add some noise to the generated fictitious odor response
            getPNStdevBhandawat(PNactivity) .* ... %returns hard-coded glomeruli-specific std-value from Bhandawat2007
            noise.*randn(nPNs, K_odors); %noise amplitude times normal distr. number
    end
        
    % making sure that PN responses are +ve
    PNtrials(PNtrials<0)=0;
    
    %rescaling the PN responses in range from 0 to 5
	%  inherited from Abdelrahman et al 2021
    PNtrials = rescale(PNtrials,0,5);

    % tune the models on odors 1:110 only - we'll use the others in the test
    % but we don't want to force average 10% coding level on a set of odors
    % that includes very extreme dilutions
    PNtrials_pure = PNtrials(:,1:110,:);
    PNtrials_pure = reshape(PNtrials_pure, size(PNtrials_pure,1),[]); %make 2-D
    
	%% begin optimisation procedure to get C_theta and APLgain values that satisfy coding level constraints
	% for the variable weight, variable theta model.

    %optimises and returns an updated MBmodel

    % random
	MBmodel_varN_varW_varTh_noHomeo = optimiseMBparams_APLgain_Ctheta(MBmodel_varN_varW_varTh, PNtrials_pure);

    % homogeneous
	% MBmodel_fixN_fixW_fixTh = optimiseMBparams_APLgain_Ctheta(MBmodel_fixN_fixW_fixTh, PNtrials_pure);

    % P(W|N,theta)
	MBmodel_varN_varW_varTh_actIndepHomeo = optimiseMBparams_APLgain_Ctheta(MBmodel_varN_varW_varTh_actIndepHomeo, PNtrials_pure);
    MBmodel_varN_varW_varTh_actIndepHomeo.modelType = 'P(w|N,theta) aka cyan model';

    % activity-dependent compensation
    MBmodel_varN_varW_varTh_withHomeo = optimiseMBparams_homeostaticExcitation(MBmodel_varN_varW_varTh, PNtrials_pure);
 
    %% save the parameters for each fly, connectivity weights, spiking thresholds, inhibition Gains,for all 8 models
    
    % save the flies tuned on original Hallem-Carlson Data:        
    allTestedOdors{randomTrials} = PNtrials; % we save the whole 176 odors so they can be used later

    % MBmodels_all{1,randomTrials} = MBmodel_fixN_fixW_fixTh;
    MBmodels_all{1,randomTrials} = MBmodel_varN_varW_varTh_noHomeo; % 'random' model
    MBmodels_all{2,randomTrials} = MBmodel_varN_varW_varTh_actIndepHomeo; % activity-independent equalization
    MBmodels_all{3,randomTrials} = MBmodel_varN_varW_varTh_withHomeo; % activity-dependent equalization
    MBmodels_all{4,randomTrials} = optimiseMBparams_homeostaticExcitation(MBmodel_varN_varW_varTh, PNtrials_pure, 0.2); % restricted activity-dependent equalization

end

save(strcat('MBmodels-',char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')),'.mat'),'MBmodels_all','allTestedOdors')

