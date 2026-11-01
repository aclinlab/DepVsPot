load('naiveParticles.mat');
out = naiveParticles;

load('fractions.mat')


[numDilutions] = length(out);
N = 2000;

numModels = 2;

modelMakeup = zeros(numModels, numDilutions); % models
meanThetas = zeros(4,numModels,numDilutions); % theta, models
hdi95s = zeros(4,2,numModels,numDilutions); % theta, lower/upper, models
for d=1:numDilutions
    for m=1:numModels
        m_indices = (out(d).theta(:,5)==m);
        modelMakeup(m,d) = sum(m_indices)/N;
        meanThetas(:,m,d) = out(d).weights(m_indices)' * out(d).theta(m_indices,1:4);
        for param=1:4
            hdi95s(param,:,m,d) = hdi95_weighted(out(d).theta(m_indices,param), out(d).weights(m_indices));
        end
    end
end

[~,winningModel] = max(modelMakeup,[],1);
winningModel = squeeze(winningModel);

bandwidthDivFactor = 2;
weighted_means = zeros(7,4);
densities = cell(4,7);
for j=1:7
    for i=1:4
        indices = out(j).theta(:,5)==winningModel(j);
        weighted_means(j,i) = mean( out(j).weights' * out(j).theta(indices,i) );
        [~,~,u] = ksdensity(out(j).theta(indices,i), 'Weights', out(j).weights);
        [densities{i,j}(:,2), densities{i,j }(:,1)] = ksdensity(out(j).theta(indices,i), 'Weights', out(j).weights,'Bandwidth',u/bandwidthDivFactor);
    end
end

muXenter_0 = mean(weighted_means(:,1));
sigmaXenter_0 = mean(weighted_means(:,2));

muXstay_shifted = weighted_means(:,3) - (-muXenter_0); % -muXenter_0 is the y-intercept
slope_muXstay_0 = fractions \ muXstay_shifted;

sigmaXstay_shifted = weighted_means(:,4) - (sigmaXenter_0); % sigmaXenter_0 is the y-intercept
slope_sigmaXstay_0 = fractions \ sigmaXstay_shifted;

%% Fig 7E

figure
ylabels = {'\mu_{enter}','\sigma_{enter}','\mu_{stay}','\sigma_{stay}'};
for i=1:4
    subplot(1,4,i), hold on

    distributionPlot(squeeze(densities(i,:)),'xValues',fractions,'showMM',0,'globalNorm',1)
    if i==1 || i ==3
        yl = [-2 2];
        yticks(-2:1:2)
        xl=xlim;
        line(xl,[0 0],'Color','k');
    else
        yl = [0 1.6];
    end
    scatter(fractions,weighted_means(:,i),2, 'MarkerFaceColor','r','MarkerEdgeColor','none')
    ylim(yl);
    xticks([0 0.1 0.2])
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    ylabel(ylabels{i})
    xlabel('Fraction KCs responding')
    axis square
    set(gca,'units','centimeters','position',[1+3*(i-1),1,2,2])
    fontsize(6,'points')

    switch i
        case 1
            plot(fractions, muXenter_0*ones(numDilutions,1),'r')
        case 2
            plot(fractions, sigmaXenter_0*ones(numDilutions,1),'r')
        case 3
            plot(fractions, -muXenter_0 + slope_muXstay_0*fractions,'r');
        case 4
            plot(fractions, sigmaXenter_0 + slope_sigmaXstay_0*fractions,'r')
    end

end

%% compare all simulated naive ecdfs vs real data

load('integV.mat')
load('naiveNullFits.mat'); % null_fits

M = 200;

simDataRealAll = cell(numDilutions,1);
simDataNullAll = cell(numDilutions,1);
for d=1:numDilutions
    theta = [weighted_means(d,:), winningModel(d)];
    theta_null = null_fits(d).weights' * null_fits(d).theta;
    theta_null([2 4]) = 0;
    decsCounts = decisionsTotalIntegV{d}(:,2);
    simDataReal = zeros(length(decsCounts),M,3);
    simDataNull = zeros(length(decsCounts),M,3);
    parfor i=1:M
        simDataReal(:,i,:) = simFlyDecs_MS(theta, decsCounts);
        simDataNull(:,i,:) = simFlyDecs(theta_null, decsCounts); % default model
    end
    simDataRealAll{d} = simDataReal;
    simDataNullAll{d} = simDataNull;
end


%% Fig S3 (7C is part of this)

simDataToPlot = [simDataNullAll,simDataRealAll];
forPub = true;
figure('Position',[100 100 600 600])
obsData = cell(numDilutions,1);
for d=1:numDilutions
    decsCounts = decisionsTotalIntegV{d}(:,2);
    numFlies = length(decsCounts);
    obsData{d} = zeros(numFlies,3);
    obsData{d}(:,1) = trainedOdorDecsIntegV{d}(:,2);
    obsData{d}(:,2) = decsTrainedOdorStayIntegV{d}(:,2);
    obsData{d}(:,3) = decsTrainedOdorEnterIntegV{d}(:,2);
    for m=1:3 %app, stay, enter

        for i=1:2 %null, real
            if forPub
                ax = axes;
            else
                subplot(numDilutions, 6, (d-1)*6 + (m-1)*2 + i);
            end
            compareECDFs(obsData{d}(:,m), simDataToPlot{d,i}(:,:,m))

            if forPub
                set(gca,'TickLength',[.03 .03],'TickDir','out')
                axis square
                set(gca,'units','centimeters','position',[3*(2*(m-1)+0.9*i),1+2.7*(7-d),2,2])
                fontsize(6,'points')
                set(gca,'color','none')
            end
        end
    end

end
