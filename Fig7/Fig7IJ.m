load('JulianMelissaData.mat') 
load('fractions.mat')

numDilutions = size(decisionsTotal,1);
numVoltages = size(decisionsTotal,2);

% Extract Observed Decision Counts for the Learned Phase

e = [2 6]; % pre post
prePostDecsCounts = cell(numDilutions,numVoltages);
for d=1:numDilutions
    for v=1:numVoltages
        prePostDecsCounts{d,v} = decisionsTotal{d,v}(:,e);
    end
end

obsData = cell(numDilutions, numVoltages);
e = [2 6];
for d = 1:numDilutions
    for v = 1:numVoltages
        numFlies = length(decisionsTotal{d,v});
        obsData{d,v} = zeros(numFlies,3,2);
        obsData{d,v}(:,1,:) = trainedOdorDecs{d,v}(:,e);
        obsData{d,v}(:,2,:) = decsTrainedOdorStay{d,v}(:,e);
        obsData{d,v}(:,3,:) = decsTrainedOdorEnter{d,v}(:,e);
    end
end


% Tm=0 particle
S = load('learnedTmZero.mat'); 
negCtrl_thetas_all = [];
negCtrl_weights_all = [];
for i=1:length(S.results)
    negCtrl_thetas_all = [negCtrl_thetas_all; S.results(i).theta];
    negCtrl_weights_all = [negCtrl_weights_all; S.results(i).weights];
end
negCtrl_weights_all = negCtrl_weights_all/sum(negCtrl_weights_all);
negCtrl = negCtrl_weights_all' * negCtrl_thetas_all;
negCtrl(14) = 0;
null_map_particle = negCtrl;

% this comes from the script FigS4E_S5.m
real_map_particle = map_particles(1,:);


%% run simulations on the particles from the null and real fit

M = 200;

simDataReal = cell(M,numDilutions,numVoltages);
simDataNull = cell(M,numDilutions,numVoltages);

parfor i=1:M
    simDataNull(i,:,:) = simFlyDecsLearnedTmK(null_map_particle,prePostDecsCounts, fractions);
    simDataReal(i,:,:) = simFlyDecsLearnedTmK(real_map_particle,prePostDecsCounts, fractions);
end
% need to combine after the parfor loop
simData = cat(4,simDataNull, simDataReal); % dims: M, numDilutions, numVoltages, 2

% each cell in simData is a matrix with numFlies x 3 x 2
% get the average across flies
sim_mean_across_flies = cellfun(@(x) squeeze(mean(x,1,'omitnan')), simData, 'UniformOutput',false);
obs_mean_across_flies = cellfun(@(x) squeeze(mean(x,1,'omitnan')), obsData, 'UniformOutput',false);

% calculate means and CI95s across repetitions
meansToPlot = zeros(numDilutions, numVoltages, 2, 3); % num dilutions, num voltages, null vs real, which metric
ci95ToPlot = zeros(numDilutions, numVoltages, 2, 2, 3); % num dilutions, num voltages, lower/upper, null vs real, which metric
targets = zeros(numDilutions, numVoltages, 3);
for d=1:numDilutions
    for v=1:numVoltages
        for m=1:3 % which metric
            targets(d,v,m) = cellfun(@(x) x(m,2), obs_mean_across_flies(d,v));
            for c=1:2 % null vs real
                % extract the values across all repetitions, for this
                % dilution, voltage, condition and metric
                values = cellfun(@(x) x(m,2), sim_mean_across_flies(:,d,v,c));
                meansToPlot(d,v,c,m) = mean(values,'omitnan');
                ci95ToPlot(d,v,c,1,m) = prctile(values, 2.5);
                ci95ToPlot(d,v,c,2,m) = prctile(values, 97.5);
            end
        end
    end
end

%% plot CvM divergence for null vs real fits (for the 200 simulations)

% load each ABC run separately

nDims = 18;
nRuns = 10;
[real_particles,null_particles] = deal(zeros(nDims,nRuns));
S = load('learnedTmZero.mat');
for i=1:length(S.results)
    null_particles(:,i) = S.results(i).weights' * S.results(i).theta;
end
S = load('learnedRealData.mat');
for i=1:length(S.results)
    real_particles(:,i) = S.results(i).weights' * S.results(i).theta;
end

% run simulations for every particle
M = 200;

simDataRealRuns = cell(nRuns,M,numDilutions,numVoltages);
simDataNullRuns = cell(nRuns,M,numDilutions,numVoltages);

for j=1:nRuns
    j
    parfor i=1:M
        simDataNullRuns(j,i,:,:) = simFlyDecsLearnedTmK(null_particles(:,j)',prePostDecsCounts, fractions);
        simDataRealRuns(j,i,:,:) = simFlyDecsLearnedTmK(real_particles(:,j)',prePostDecsCounts, fractions);
    end
end
% need to combine after the parfor loop

simDataRuns = cat(5,simDataNullRuns, simDataRealRuns); % dims: M, numDilutions, numVoltages, 2
% 1 is null
% 2 is real

lastEpoch = cellfun(@(x) x(:,:,2), simDataRuns, 'UniformOutput', false); % slice the last epoch

%% calc CvM divergence

[interSimCMD, dataVsSimCMD] = ...
    deal(nan(numDilutions,numVoltages,3,2,nRuns)); % 3 = num measures; 2 = real vs null
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s\n', msg));
parfor j=1:nRuns
    for d=2:numDilutions
        for v=2:numVoltages
            send(q, sprintf('run %d dilution %d voltage %d\n',j,d,v));
            for c=1:2 % real vs null
                [~, ~, dataVsSim, interSim] = ...
                    MonteCarloCM2test(obsData{d,v}(:,:,2), permute(cat(3, lastEpoch{j,:,d,v,c}), [1 3 2]));
                dataVsSimCMD(d,v,:,c,j) = dataVsSim;
                interSimCMD(d,v,:,c,j) = mean(interSim,1);
            end
        end
    end
end

%% Fig 7I: plot CvM divergence
dataVsSim = squeeze(mean(dataVsSimCMD(2:end,2:end,:,:,:),[1 2],'omitnan')); % measures, null vs real, particles
interSim = squeeze(mean(interSimCMD(2:end,2:end,:,:,:),[1 2],'omitnan')); % measures, null vs real, particles
for m=[1 3]
    groupNames  = {'T_m = 0', 'Fitted T_m'};
    seriesData  = {squeeze(dataVsSim(m,:,:)), squeeze(interSim(m,:,:))};
    seriesNames = {'data vs sim', 'inter-sim'};
    seriesCols  = {[0 0.4470 0.7410], [0 0 0]};   % blue [default first color for Matlab], black
    markerSize  = 10;                       % scatter area (points^2)
    jitterWidth = 0.25;                     % total horizontal spread of each swarm
    medHalfW    = 0.12;                     % half-width of median line

    nGroups = size(dataVsSim, 2);

    figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 10 6.5]);
    ax = axes; hold(ax, 'on');

    hLeg = gobjects(1, numel(seriesData));
    for s = 1:numel(seriesData)
        Y = seriesData{s};
        for g = 1:nGroups
            y = Y(g, :);
            x = g * ones(size(y));
            hs = swarmchart(ax, x, y, markerSize, seriesCols{s}, 'filled', ...
                'XJitter', 'density', 'XJitterWidth', jitterWidth,...
                'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
            if g == 1, hLeg(s) = hs; end

            meany = mean(y, 'omitnan');
            plot(ax, g + [-1 1]*medHalfW, [meany meany], '-', ...
                'Color', seriesCols{s}, 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end

    % --- axis styling ---
    set(ax, 'XTick', 1:nGroups, 'XTickLabel', groupNames, ...
        'XLim', [0.5 nGroups + 0.5], ...
        'TickDir', 'out', 'Box', 'off', 'LineWidth', 0.5, ...
        'FontName', 'Arial', 'FontSize', 6, 'TickLength',[.03 .03],...
        'units','centimeters','position',[1 1 2 2]);
    ax.XAxis.TickLength = [0 0];     % Prism hides x ticks for categories

    % legend(hLeg, seriesNames, 'Location', 'northeastoutside', ...
    %     'Box', 'off', 'FontWeight', 'normal');
end


%% Fig S4F: plot data vs sim for P(enter)post

plotModelComparisons(targets(:,:,3), meansToPlot(:,:,:,3), ci95ToPlot(:,:,:,:,m),'P(enter)post',true);

%% Fig 7J: plot HDI95 of posterior theta - single plot
figure
    hold on
d = 14;
conditions = [2 1];
colors = {'r','k'};
alpha = 0.2;
transparent_colors = [...
    1 0 0 alpha; % r
    0 0 0 alpha]; % k
for ci=1:length(conditions)
    c = conditions(ci);
    numReps = size(thetas_all{c},1);
    for i=1:numReps
        [f,xi] = ksdensity(thetas_all{c}(i,:,d),'Weights',weights_all{c}(i,:),'BoundaryCorrection','reflection');
        plot(xi,f,'Color',transparent_colors(c,:),'LineWidth',0.25)
    end
    [f,xi] = ksdensity(reshape(thetas_all{c}(:,:,d),[],1),'Weights',weights_all{c}(:),'BoundaryCorrection','reflection');
    plot(xi,f,colors{c},'LineWidth',1)
    yl = ylim;
    errorbar(map_particles(c,d), yl(2)-(range(yl))/10, ...
        map_particles(c,d) - hdi95_combined{c}(d,1), hdi95_combined{c}(d,2) - map_particles(c,d), ...
        'horizontal', 'o','Color',colors{c},'MarkerFaceColor',colors{c},'MarkerSize',3)
    line([0 0], yl, 'Color', 'k')
    xlabel('T_m');
    ylabel('Posterior probability density')
    xlim([-0.2 0.4])


end
ylim(yl)
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(gca,6,'points')


%% Fig S8A compare pre vs post correlations

trueCorrelations = zeros(numDilutions,numVoltages,3);
simCorrelations = zeros(numDilutions,numVoltages,3,M);
corrZscores = zeros(numDilutions,numVoltages,3);
avgCorrDiff = zeros(M,3);
normAvgCorrDiff = zeros(M,1);
for d=1:numDilutions
    for v=1:numVoltages
        for w = 1:3 % which metric (approach, stay, enter)
            nans = isnan(obsData{d,v}(:,w,1)) | isnan(obsData{d,v}(:,w,2));
            trueCorrelations(d,v,w) = corr(obsData{d,v}(~nans,w,1),obsData{d,v}(~nans,w,2),'Type','Kendall');% corr(squeeze(obsData{d,v}(:,1,:)));
            for i=1:M
                nans = isnan(simDataReal{i,d,v}(:,w,1)) | isnan(simDataReal{i,d,v}(:,w,2));
                simCorrelations(d,v,w,i) = corr(simDataReal{i,d,v}(~nans,w,1),simDataReal{i,d,v}(~nans,w,2),'Type','Kendall');%corr(squeeze(simDataReal{i,d,v}(:,1,:)));
            end
            corrZscores(d,v,w) = (trueCorrelations(d,v,w) - mean(simCorrelations(d,v,w,:))) / std(simCorrelations(d,v,w,:));
        end
    end
end
for i=1:M
    for w=1:3
        avgCorrDiff(i,w) = mean(simCorrelations(:,:,w,i),'all') - mean(trueCorrelations(:,:,w),'all');
    end
    normAvgCorrDiff(i) = norm(avgCorrDiff(i,:));
end
forPub = true;
if forPub
    figure('Name','Kendall')
    % for w=1:3
    w=1;
        for d=1:numDilutions
            for v=1:numVoltages
                % subplot(numVoltages*3, numDilutions,d+(v-1+(w-1)*numVoltages)*numDilutions)
                if forPub
                    ax = axes;
                else
                    subplot(numDilutions,numVoltages,v+(d-1)*numVoltages)
                end
                hold on
                % histogram(squeeze(simCorrelations(d,v,:)))
                [f,xi] = ksdensity(squeeze(simCorrelations(d,v,w,:)));
                plot(xi,f)
                line([0 0], [0 10],'Color','k')
                line(trueCorrelations(d,v,w)*[1 1], [0 10],'Color','r')
                xlim([-0.2 0.6])
                ylim([0 10])
                axis square
                set(gca,'TickLength',[.03 .03],'TickDir','out')
                if d~=7
                    set(gca,'XTickLabel',[])
                end
                if v~=1
                    set(gca,'YTickLabel',[])
                end
                set(gca,'units','centimeters','position',[1.2*v+1,1+1.2*(7-d),1,1])
                fontsize(6,'points')
                set(gca,'color','none')

            end
        end
    % end
else
    for w=1:3
        figure('Name',['Kendall ' num2str(w)])
        for d=1:numDilutions
            for v=1:numVoltages
                subplot(numDilutions,numVoltages,v+(d-1)*numVoltages)
                hold on
                % histogram(squeeze(simCorrelations(d,v,:)))
                [f,xi] = ksdensity(squeeze(simCorrelations(d,v,w,:)));
                plot(xi,f)
                line([0 0], [0 10],'Color','k')
                line(trueCorrelations(d,v,w)*[1 1], [0 10])
                xlim([-0.2 0.6])
                ylim([0 10])
            end
        end
    end
end

%% Fig S8B example correlations
figure,scatter(obsData{5,1}(:,1,1),obsData{5,1}(:,1,2),5)
xlabel('p_{obs}(approach)_{pre}')
ylabel('p_{obs}(approach)_{post}')
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(gca,6,'points')

%% Fig S8c: plot ks density of k


figure
d = 18;
ci = 2;
conditions = [2 1]; % just plot the real data
% for ci=1:length(conditions)
    c = conditions(ci);
    axHandles(ci) = subplot(length(conditions),1,ci);
    hold on
    numReps = size(thetas_all{c},1);
    for i=1:numReps
        [f,xi] = ksdensity(1-thetas_all{c}(i,:,d),'Weights',weights_all{c}(i,:),'BoundaryCorrection','reflection');
        plot(xi,f,'Color',transparent_colors(c,:),'LineWidth',0.25)
    end
    [f,xi] = ksdensity(1-reshape(thetas_all{c}(:,:,d),[],1),'Weights',weights_all{c}(:),'BoundaryCorrection','reflection');
    plot(xi,f,colors{c},'LineWidth',1)
    yl = ylim;
    errorbar(1-map_particles(c,d), yl(2)-(range(yl))/10, ...
        1-map_particles(c,d) - (1-hdi95_combined{c}(d,1)), (1-hdi95_combined{c}(d,2)) - (1-map_particles(c,d)), ...
        'horizontal', 'o','Color',colors{c},'MarkerFaceColor',colors{c},'MarkerSize',3)
    xlabel('k');
    ylabel('P(k|data) density')
    xlim([0 1])

    set(axHandles(ci),'units','centimeters','position',[1,7-3*(ci-1),2,2])
    set(axHandles(ci),'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(axHandles(ci),6,'points')

% end



%% Fig S7: plot ecdfs simulation vs data


% simDataToPlot = [simDataNullAll,simDataRealAll];
forPub = true;
figure('Position',[100 100 600 500])
% tiledlayout(7,12);
% tiledlayout(12,7);
tiledlayout(7,8);
for d=1:numDilutions
    for m=3 %app, stay, enter

        for v = 1:numVoltages
            numFlies = size(obsData{d,v},1);
            thisSimDataReal = zeros(numFlies,M);
            thisSimDataNull = zeros(numFlies,M);
            for i=1:M % num simulations
                thisSimDataReal(:,i) = simDataReal{i,d,v}(:,m,2); % post data
                thisSimDataNull(:,i) = simDataNull{i,d,v}(:,m,2); % post data
            end
            nexttile((d-1)*(numVoltages*2) + v)
            compareECDFs(obsData{d,v}(:,m,2), thisSimDataNull)

            if forPub
                set(gca,'TickLength',[.03 .03],'TickDir','out')
                axis square
                % set(gca,'units','centimeters','position',[1+2.7*(d-1),3*(2*()+0.9*m),2,2])
                fontsize(6,'points')
                % set(gca,'color','none')
            end

            nexttile((d-1)*(numVoltages*2) + v + numVoltages)
            compareECDFs(obsData{d,v}(:,m,2), thisSimDataReal)
            if forPub
                set(gca,'TickLength',[.03 .03],'TickDir','out')
                axis square
                % set(gca,'units','centimeters','position',[1+2.7*(d-1),3*(2*()+0.9*m),2,2])
                fontsize(6,'points')
                % set(gca,'color','none')
            end

        end
    end
end

