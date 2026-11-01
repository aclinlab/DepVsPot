 clear all

numtrainingSamples=15;

n=2000;
m=24;

logLRs = [-5 -4 -3.5 -3 -2.75 -2.5 -2.25 -2 -1.75 -1.25 -1 -0.5 0]; %-1.25; %
lrs=size(logLRs,2);

odors=176; %mulOd;

% Upper bound for the KC-MBON weights =1
Wmax=1.*ones(2000,1);

% load PNs
load('carlson_frombase.mat','carlson_frombase');
PN_frombase = ORNtoPNmodelolsen(carlson_frombase);
PN_frombase = PN_frombase'; %numPNs x numOdors
numTrials = 30;
pureOdors = 1:110;
PN_frombase_all = PN_frombase;
PN_frombase_pure = PN_frombase(:,pureOdors);

% The order of models in the matrix goes like:
% 1. Homogeneous (black) - all parameters equal - formerly no suffix
% 2. Random (red) - all variable, no compensation - formerly "Homog
% Fixedtheta)
% 3. Activity-independent compensation (cyan) - "Equalized"
% 4. Adjust weights (blue) - "comp2_noxjk"
% This is the order they were presented in, in Abdelrahman et al 2021 PNAS

load('MBmodels.mat','MBmodels_all','allTestedOdors');

numModelInstances=size(MBmodels_all,2);
numModels = size(MBmodels_all,1); 
WmaxExp = repmat(Wmax,1,numModels,odors);

test_p = zeros(numModelInstances, lrs, numModels, odors, odors, 2);
test_p_gen = zeros(numModelInstances, lrs, numModels, odors, odors, 2);

Ys_binarised_all = false(numModelInstances,n,odors,numTrials,numModels);

% if only running one model instance, don't use parfor
parfor randomTrials=1:numModelInstances
    randomTrials

    Ys = zeros(n,odors,numTrials,numModels);
    for mi = 1:numModels
        Ys(:,:,:,mi) = reshape(calculateKCresponse(MBmodels_all{mi,randomTrials},reshape(allTestedOdors{randomTrials},m,[])),n,odors,numTrials);
    end

    % binarise KC activity
    Ys_binarised = (Ys>0); 

    Ytrain = Ys_binarised(:,:,1:numtrainingSamples,:);

    this_test_p = zeros(1,lrs,numModels,odors,odors,2);
    this_test_p_gen = zeros(1,lrs,numModels,odors,odors,2);
    for l_r=1:lrs
        l_r
        tic

        % Wop will be: numKCs x 2 x numModels
        Wop = repmat(zeros(n,2),1,1,numModels, odors, 2); %last dim is for dep vs pot
        Wop(:,:,:,:,2) = 1; % for depression, start at 1

        alpha=10.^logLRs(l_r);
        alphaNormToModel = alpha./squeeze(mean(Ytrain,[1 2 3])); % avg across all KCs, odors and trials to get avg Y for each model
        alphaNormToModel = repmat(alphaNormToModel',n,1); % will end up as numKCs x numModels

        %% Training a preceptron in the output layer: learning by potentiating the right action

        % this is a vectorized implementation
        % you can sum up all odors of the same valence and combine 

        % for the odor x odor square
        alphaNormToModel = repmat(alphaNormToModel, 1, 1, odors);
        alphaNormToModel = permute(alphaNormToModel, [1 3 2]); % dimensions are now numKCs x numOdors x numModels
        deltaByOdor = exp( -alphaNormToModel .* squeeze( sum(Ytrain, 3) ) ); % numKCs x numOdors x numModels
        deltaByOdor = permute(deltaByOdor, [1 3 2]); % numKCs x numModels x numOdors

        % potentiation
        % the reason for using "permute" is to squeeze out only the 2nd
        % dimension
        Wop(:,1,:,:,1) = WmaxExp - (WmaxExp-permute(Wop(:,1,:,:,1),[1 3 4 5 2])).*deltaByOdor;
        Wop(:,2,:,:,1) = WmaxExp - (WmaxExp-permute(Wop(:,2,:,:,1),[1 3 4 5 2])).*deltaByOdor;

        % depression
        Wop(:,1,:,:,2) = permute(Wop(:,1,:,:,2),[1 3 4 5 2]).*deltaByOdor;
        Wop(:,2,:,:,2) = permute(Wop(:,2,:,:,2),[1 3 4 5 2]).*deltaByOdor;

        for model = 1:numModels
            model
            for odor1 = 1:odors
                for odor2 = 1:odors

                    % discrimination 

                    % potentiation
                    [this_test_p(1, l_r, model, odor1, odor2,1),~] = ... 
                        testing_models_accuracies_normcdf([Wop(:,1,model,odor1,1), Wop(:,2,model,odor2,1)], ...
                        1, numTrials, numtrainingSamples, Ys_binarised(:,odor1,:,model), 'div');

                    % depression
                    [this_test_p(1, l_r, model, odor1, odor2,2),~] = ... 
                        testing_models_accuracies_normcdf([Wop(:,1,model,odor1,2), Wop(:,2,model,odor2,2)], ...
                        2, numTrials, numtrainingSamples, Ys_binarised(:,odor1,:,model), 'div');

                    % generalization

                    % potentiation
                    [this_test_p_gen(1, l_r, model, odor1, odor2,1),~] = ... 
                        testing_models_accuracies_normcdf([Wop(:,1,model,odor1,1), zeros(n,1)], ...
                        1, numTrials, numtrainingSamples, Ys_binarised(:,odor2,:,model), 'div');
                    % depression
                    [this_test_p_gen(1, l_r, model, odor1, odor2,2),~] = ... 
                        testing_models_accuracies_normcdf([ones(n,1), Wop(:,1,model,odor1,2)], ...
                        1, numTrials, numtrainingSamples, Ys_binarised(:,odor2,:,model), 'div');
                    
                end
            end
        end
        toc

        test_p(randomTrials,:,:,:,:,:) = this_test_p;
        test_p_gen(randomTrials,:,:,:,:,:) = this_test_p_gen;
        Ys_binarised_all(randomTrials,:,:,:,:) = Ys_binarised;

    end
end

% convert P(enter) to P(reverse)
test_p = 1 - test_p;
test_p_gen = 1 - test_p_gen;

%%
save(['LearningByDepPotPairs_results',char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')),'.mat'],...
    'test_p','test_p_gen','logLRs','Ys_binarised_all', 'odors', 'numModels');
