%% load data
% (not in repo because of file sizes but equivalent files will be generated
% by Fig7H_unifiedModel.m, replace these with new file names)
fileList = {'unifiedModelResults260716T234824.mat';
'unifiedModelResults260717T071332.mat';
'unifiedModelResults260717T064851.mat';
'unifiedModelResults260717T062331.mat';
'unifiedModelResults260717T055837.mat';
'unifiedModelResults260717T053339.mat';
'unifiedModelResults260717T050919.mat';
'unifiedModelResults260717T044417.mat';
'unifiedModelResults260717T041940.mat';
'unifiedModelResults260717T035437.mat';
'unifiedModelResults260717T033035.mat';
'unifiedModelResults260717T030545.mat';
'unifiedModelResults260717T024056.mat';
'unifiedModelResults260717T021616.mat';
'unifiedModelResults260717T015159.mat';
'unifiedModelResults260717T012711.mat';
'unifiedModelResults260717T010217.mat';
'unifiedModelResults260717T003733.mat';
'unifiedModelResults260717T001232.mat';
'unifiedModelResults260716T232237.mat'};

Tms = [0, 0.1, 0.15, 0.175, 0.2, 0.225, 0.25, 0.3];
dataVsSimDivAll = zeros(3, length(Tms), length(fileList));
withinSimDivAll = zeros(3, length(Tms), length(fileList));
allReplicatesThetas = zeros(length(fileList),length(Tms),14);
for i=1:length(fileList)
    S = load(fileList{i});
    for j=1:length(S.paramsByTmCond)
        if S.paramsByTmCond(j).Tm ~= Tms(j)
            error('mismatched Tm')
        end
        allReplicatesThetas(i,j,5:7) = S.paramsByTmCond(j).muSs(2:4); % mu_S 15, 30, 60 V
        allReplicatesThetas(i,j,8:11) = S.paramsByTmCond(j).offsets; % H_0 at each voltage
        allReplicatesThetas(i,j,12) = S.paramsByTmCond(j).H;
        allReplicatesThetas(i,j,13) = S.paramsByTmCond(j).Ta;
        allReplicatesThetas(i,j,14) = S.paramsByTmCond(j).Tm;
        allReplicatesThetas(i,j,15:17) = S.paramsByTmCond(j).sigmaSs(2:4);
    end
    dataVsSimDivAll(:,:,i) = S.dataVsSimDiv;
    withinSimDivAll(:,:,i) = S.withinSimDiv;
end

%% figure
dataVsSimDiv_ci95 = ConfInt95(dataVsSimDivAll,3);
withinSimDiv_ci95 = ConfInt95(withinSimDivAll,3);
dataVsSimDiv_std = std(dataVsSimDivAll,0,3);
withinSimDiv_std = std(withinSimDivAll,0,3);
dataVsSimDiv_mean = mean(dataVsSimDivAll,3);
withinSimDiv_mean = mean(withinSimDivAll,3);
figure
subplotLabels = {'P(appr)','P(stay)','P(enter)'};
legendLabels = {'data vs. sim', 'sim vs. sim'};
for p=[1 3] % figure plots P(appr) and P(enter)
    subplot(1,3,p)
    title(subplotLabels{p})
    hold on
    errorbar(Tms,dataVsSimDiv_mean(p,:),dataVsSimDiv_std(p,:),'o-','MarkerSize',3,'CapSize',3)
    errorbar(Tms,withinSimDiv_mean(p,:),withinSimDiv_std(p,:),'o-','MarkerSize',3,'CapSize',3)
    xlabel('T_m')
    ylabel('Average C-vM divergence')
    % legend(legendLabels)
    set(gca,'units','centimeters','position',[1+3*(p-1),1,2,2])
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(6,'points')
    xticks(0:0.1:0.3)
end
% totalDiffDiv = sum(dataVsSimDiv-withinSimDiv,1);
% [~, bestTmIndex] = min(totalDiffDiv);

%% plot theta values
% The range of parameter values resulting from gradient descent were used
% to guide the priors for these parameters for ABC
% this section plots the parametres values from gradient descent
% within the range eventually used for the prior
dimensionsToPlot = [5:13, 15:17];
bayesian_priors = [
    -1 0; % muS 15 
    -2 0; % muS 30 
    -3 0; % muS 60 
    -0.3 0.3; % offset 0
    -0.3 0.3; % offset 15
    0 0.6; % offset 30
    0 0.7; % offset 60
    -30 -6;  % H
    0    0.3; %Ta
    -0.5 0.5; %Tm
    0 1.5; % sigmaS 15 
    0 2; % sigmaS 30 
    0 3; % sigmaS 60 
];
paramlabels = {'\mu_{enter,pre}', '\sigma_{enter,pre}', 'm_{\mu,stay,pre}'...
    'm_{\sigma,stay,pre}','\mu_s(V=15)', '\mu_s(V=30)', '\mu_s(V=60)',...
    'H_0(V=0)','H_0(V=15)','H_0(V=30)','H_0(V=60)', 'H',...
    'T_a','T_m','\sigma_s(V=15)','\sigma_s(V=30)','\sigma_s(V=60)',...
    'unstable variability'};

figure
for i=1:length(dimensionsToPlot)
    subplot(3,4,i);
    errorbar(Tms,mean(allReplicatesThetas(:,:,dimensionsToPlot(i)), 1), ...
        std(allReplicatesThetas(:,:,dimensionsToPlot(i)), 0, 1),'o-','MarkerSize',3,'CapSize',3)
    xlabel('T_m')
    ylabel(paramlabels{dimensionsToPlot(i)})
    line([0 0.25], bayesian_priors(dimensionsToPlot(i)-4,1),'Color','k')
    line([0 0.25], bayesian_priors(dimensionsToPlot(i)-4,2),'Color','k')
end