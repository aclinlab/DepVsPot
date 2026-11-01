%% concatenate data files

% Data files are from Hige et al 2015 Neuron
% files{1} = neg controls
%   1I: US comes before CS
%   S5C: no CsChrimson in DAN
% files{2} = learning
%   1E: OCT CS+, MCH CS-
%   4D: OCT CS+, MCH CS- (with voltage clamp during training)
%   5G: OCT CS+, MCH CS- with TH-GAL4 instead of PPL1-g1-specific driver
%   S6D: OCT CS+, MCH CS- with 1 min odor+DAN pairing instead of 1 sec
%   S6H: OCT CS+, MCH CS- with only one light pulses instead of 4
%   S3C: MCH CS+, OCT CS-
%   7D: PA CS+ (also causes depression in similar odors BA,HP), EL CS-
%   7G: EL CS+, PA/BA/HP CS-
% data is organised as:
% data{condition}(cells) - where each cell is a struct containing:
%  CSplus
%  CSminus
% learned data only:
%  CSp1
%  CSp2
%  CSp3
%  CSm1
%  CSm2
%  CSm3


files = { {'Fig1I.mat', 'FigS5C.mat'}, ...
    {'Fig1E.mat',
    'Fig4D.mat',
    'Fig5G.mat',
    'FigS6D.mat',
    'FigS6H.mat',
    'FigS3C.mat',
    'Fig7D.mat',
    'Fig7G.mat'} };

odorNames = {'CSplus','CSminus','CSp1','CSp2','CSp3','CSm1','CSm2','CSm3'};
numOdors = length(odorNames);

forFigure = true; % for formatting the figures into publication style

clear data
for condition=1:2 % control vs learning
    count = 1;
    for f=1:length(files{condition})
        thisData = load(files{condition}{f});
        fieldNames = fieldnames(thisData);
        cellsInThisFile = length(thisData.(fieldNames{1})); % assume all odors have the same number of cells

        % Combine all conditions regardless of which odor was CS+/CS-
        switch files{condition}{f}
            case {'Fig1E.mat','Fig4D.mat','Fig5G.mat','FigS6D.mat','FigS6H.mat','Fig1I.mat','FigS5C.mat'}
                for i=1:cellsInThisFile
                    data{condition}(count+i-1).CSplus = thisData.OCT{i};
                    data{condition}(count+i-1).CSminus = thisData.MCH{i};
                end
            case 'FigS3C.mat'
                for i=1:cellsInThisFile
                    data{condition}(count+i-1).CSplus = thisData.MCH{i};
                    data{condition}(count+i-1).CSminus = thisData.OCT{i};
                end
            case 'Fig7D.mat'
                for i=1:cellsInThisFile
                    data{condition}(count+i-1).CSp1 = thisData.PA{i};
                    data{condition}(count+i-1).CSminus = thisData.EL{i};
                    data{condition}(count+i-1).CSp2 = thisData.BA{i};
                    data{condition}(count+i-1).CSp3 = thisData.HP{i};
                end
            case 'Fig7G.mat'
                for i=1:cellsInThisFile
                    data{condition}(count+i-1).CSplus = thisData.EL{i};
                    data{condition}(count+i-1).CSm1 = thisData.PA{i};
                    data{condition}(count+i-1).CSm2 = thisData.BA{i};
                    data{condition}(count+i-1).CSm3 = thisData.HP{i};
                end
        end
        count = count + cellsInThisFile;

    end
end

%% set parameters

% Parameters
odor_onset = 3.2;
bin_size = 0.05; % 50 ms bins for spike rate
time_window = [0 10]; % take the full 10 s
edges = time_window(1):bin_size:time_window(2);
max_time = 10;
numBins = length(edges)-1;
binCenters = edges(1:(end-1))+bin_size/2;

% split the recordings into 3 periods:
% 0 - 3.35 s: pre odor - although odor onset is at 3.2 s, it's visually
% obvious from inspecting the spike trains that the response doesn't start until 3.35 s
% 3.35 - 4.6 s: during or just after the odor
% 4.6 - 10 s: post odor
periodTimeEdges = [0 3.35 4.6 10];
numPeriods = length(periodTimeEdges)-1;
periods = zeros(numPeriods,2);
for p=1:numPeriods
    periods(p,1) = round(periodTimeEdges(p)/bin_size) + 1;
    periods(p,2) = round(periodTimeEdges(p+1)/bin_size);
end
numPeriods = size(periods,1);

%% Fig 4C

odorData = {data{2}.CSplus};

figure
cellIndex = 6; % the particular cell illustrated in Fig 4C
% pre/post
trialCount = 0;
for j=1:2
    % cells
    for t=1:5
        trialCount = trialCount + 1;
        spikeTimes = odorData{cellIndex}{j}{t};
        numSpikes = length(spikeTimes);
        for spike = 1:numSpikes
            line([spikeTimes(spike), spikeTimes(spike)], ...
                -[trialCount-0.4, trialCount+0.4],'Color','k');
        end
    end
    line([0 10], -[trialCount+0.5, trialCount+0.5]);
end
xlim([3 5])
exportgraphics(gcf, 'spikes.pdf', 'ContentType', 'vector')


%% Get the spike counts and inter-trial variance in each bin across all cells

% 

numConds = 2;
[variances, means, variances_by_period, means_by_period, all_trial_counts, counts_by_period] = deal(cell(numConds,1));
for c = 1:numConds
    odorNames = fieldnames(data{c}(1));
    numCells = length(data{c}); 
    numOdors = length(odorNames);
    variances{c} = nan(numOdors, 2, numBins, numCells);
    means{c} = nan(numOdors, 2, numBins, numCells);
    [variances_by_period{c}, means_by_period{c}] = deal(nan(numOdors, 2, numPeriods, numCells));
    all_trial_counts{c} = cell(numOdors, 2, numCells);
    counts_by_period{c} = cell(numOdors, 2, numCells);
    % Now calculate inter-trial variance in spike rate
    for i = 1:numCells
        for o = 1:length(odorNames)
            odorData = data{c}(i).(odorNames{o});

            if ~isempty(odorData)
                % Loop over cells
                for j = 1:2  % 1 = pre-pairing, 2 = post-pairing
                    numTrials = length(odorData{j});
                    all_trial_counts{c}{o,j,i} = zeros(numTrials,numBins);
                    counts_by_period{c}{o,j,i} = zeros(numTrials,numPeriods);

                    for t = 1:numTrials
                        % Get spike times and compute spike rate
                        spike_times = odorData{j}{t};
                        spike_times = spike_times(spike_times >= time_window(1) & spike_times <= time_window(2));

                        % Compute spike counts in bins
                        spike_counts = histcounts(spike_times, edges);
                        all_trial_counts{c}{o,j,i}(t,:) = spike_counts;
                        for p=1:numPeriods
                            counts_by_period{c}{o,j,i}(t,p) = sum(all_trial_counts{c}{o,j,i}(t,periods(p,1):periods(p,2)));% histcounts(spike_times,[periods(:,1); numBins]*bin_size);
                        end
                    end

                    % Compute trial-by-trial variance in spike rate
                    % if ~isempty(all_trial_counts)
                    var_trial_counts = var(all_trial_counts{c}{o,j,i}, 0, 1); % Variance across trials
                    mean_trial_counts = mean(all_trial_counts{c}{o,j,i}, 1); % Mean across trials
                    fano_factor = var_trial_counts ./ mean_trial_counts; % Variance-to-mean ratio

                    variances{c}(o,j,:,i) = var_trial_counts;
                    means{c}(o,j,:,i) = mean_trial_counts;

                    variances_by_period{c}(o,j,:,i) = var(counts_by_period{c}{o,j,i},0,1);
                    means_by_period{c}(o,j,:,i) = mean(counts_by_period{c}{o,j,i},1);
                    % end

                end

            end
        end

        if length(odorNames)>2
            % If it's from Fig 7, average together CSp1-3 or CSm1-3
            avgMultiCS = false;
            if ~isnan(variances{c}(contains(odorNames,'CSp1'),1,1,i)) && isnan(variances{c}(contains(odorNames,'CSplus'),1,1,i))
                target = contains(odorNames,'CSplus');
                source = contains(odorNames,{'CSp1','CSp2','CSp3'});
                avgMultiCS = true;
            elseif ~isnan(variances{c}(contains(odorNames,'CSm1'),1,1,i)) && isnan(variances{c}(contains(odorNames,'CSminus'),1,1,i))
                target = contains(odorNames,'CSminus');
                source = contains(odorNames,{'CSm1','CSm2','CSm3'});
                avgMultiCS = true;
            end
            if avgMultiCS
                variances{c}(target,:,:,i) = mean(variances{c}(source,:,:,i),1);
                means{c}(target,:,:,i) = mean(means{c}(source,:,:,i),1);
                variances_by_period{c}(target,:,:,i) = mean(variances_by_period{c}(source,:,:,i),1);
                means_by_period{c}(target,:,:,i) = mean(means_by_period{c}(source,:,:,i),1);
            end
        end
    end
    stds{c} = sqrt(variances{c});
    std_by_period{c} = sqrt(variances_by_period{c});

    FanoFactor_by_period{c} = variances_by_period{c}./means_by_period{c};
end


%% plot mean and variance across bins (Fig 4D-I)
figure

% loop over control vs. learning
for c=1:numConds
    numCells = size(means{c},4);
    meanMeanAcrossCells = mean(means{c},4);
    meanMeanAcrossCellsCI = std(means{c},0,4)/sqrt(numCells)*tinv(0.975,numCells-1);
    meanVarianceAcrossCells = mean(variances{c},4);
    meanStdAcrossCells = mean(stds{c},4);
    meanStdAcrossCellsCI = std(stds{c},0,4)/sqrt(numCells)*tinv(0.975,numCells-1);

    % loop over CS+ vs CS-
    numOdorsToPlot = 2;
    for o=1:numOdorsToPlot
        % first row (Fig 4D,F): spike rates
        subplot(3, numConds*numOdorsToPlot, o+(c-1)*2)
        hold on
        shadedErrorBar(edges(1:end-1)-3.2, squeeze(meanMeanAcrossCells(o,1,:))/bin_size,meanMeanAcrossCellsCI(o,1,:)/bin_size,'lineprops','-k');
        shadedErrorBar(edges(1:end-1)-3.2, squeeze(meanMeanAcrossCells(o,2,:))/bin_size,meanMeanAcrossCellsCI(o,2,:)/bin_size,'lineprops','-r');
        ylabel('Spikes/s')
        xlabel('Time (s)')
        xlim([-1 3])
        ylim([0 150])
        if forFigure
            set(gca,'units','centimeters','position',[1+2.8*o+(c-1)*5.9,10.5,1.52,2])
            set(gca,'TickLength',[.03 .03],'TickDir','out')
            fontsize(gca,6,'points')
        end

        % 2nd row (Fig 4E,G): std in spike rates per bin across trials
        subplot(3, numConds*numOdorsToPlot,numConds*numOdorsToPlot + o+(c-1)*2)
        hold on
        shadedErrorBar(edges(1:end-1)-3.2, squeeze(meanStdAcrossCells(o,1,:))/bin_size,meanStdAcrossCellsCI(o,1,:)/bin_size,'lineprops','-k');
        shadedErrorBar(edges(1:end-1)-3.2, squeeze(meanStdAcrossCells(o,2,:))/bin_size,meanStdAcrossCellsCI(o,2,:)/bin_size,'lineprops','-r');
        ylabel('Inter-trial st. dev. (spikes/s)')
        xlim([-1 3])
        ylim([0 30])
        xlabel('Time (s)')
        if forFigure
            set(gca,'units','centimeters','position',[1+2.8*o+(c-1)*5.9,7,1.52,2])
            set(gca,'TickLength',[.03 .03],'TickDir','out')
            fontsize(gca,6,'points')
        end

        % 3rd row (Fig H,I): std vs mean in spike count during odor
        subplot(3, numConds*numOdorsToPlot,numConds*numOdorsToPlot*2 + o+(c-1)*2)
        hold on
        for i=1:length(means_by_period{c}(o,1,2,:))
            line(squeeze(means_by_period{c}(o,:,2,i)),squeeze(std_by_period{c}(o,:,2,i)),'Color',[0 0 0 0.2])
        end
        scatter(squeeze(means_by_period{c}(o,1,2,:)), squeeze(std_by_period{c}(o,1,2,:)), 5, 'k','MarkerEdgeAlpha',0.4);%,'filled');%,'MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4)
        scatter(squeeze(means_by_period{c}(o,2,2,:)), squeeze(std_by_period{c}(o,2,2,:)), 5, 'r','MarkerEdgeAlpha',0.4);%,'filled');%,'MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4)
        xlabel('Mean spike count')
        ylabel('Inter-trial st. dev spike count')
        xlim([0 200])
        ylim([0 30])
        if forFigure
            set(gca,'units','centimeters','position',[1+2.8*o+(c-1)*5.9,3.5,1.52,1.52])
            set(gca,'TickLength',[.03 .03],'TickDir','out')
            fontsize(gca,6,'points')
        end
    end
end

%% Fig 4J-L

% get the pre vs post differences in spike count

stdevs_preVsPostDiff = cell(numConds,1);
for c=1:numConds
    stdevs_preVsPostDiff{c} = diff(std_by_period{c},1,2);
end

% std pre vs post, CS+ and CS-
if forFigure
    meanLineThickness = 1;
    circleSize = 10;
else
    meanLineThickness = 2;
    circleSize = 25;
end
figure

% Fig 4J,K
for c=1:numConds
    subplot(1,numConds+1,c)
    hold on
    meanLineWidth = 0.3;
    numCells = size(std_by_period{c},4);
    for o=1:2
        for i=1:numCells
            line([1 2]+(o-1)*2, std_by_period{c}(o,:,2,i), 'Color',[0 0 0 0.1])
            scatter([1 2]+(o-1)*2, std_by_period{c}(o,:,2,i), circleSize, 'k','MarkerEdgeAlpha',0.1)
        end
        for j=1:2
            line([j-meanLineWidth j+meanLineWidth]+(o-1)*2, ...
                [mean(std_by_period{c}(o,j,2,:),4) mean(std_by_period{c}(o,j,2,:),4)], 'Color', [0 0.5 1] , 'LineWidth',meanLineThickness)
        end
    end
    xlim([0.5 4.5])
    xticks(1:4);
    xticklabels({'Pre','Post','Pre','Post'})
    ylabel('Inter-trial st. dev. in spike count')
    ylim([0 30])
    if forFigure
        set(gca,'units','centimeters','position',[1+3.867*(c-1),1.5,2.2,2.5])
        set(gca,'TickLength',[.03 .03],'TickDir','out')
        fontsize(gca,6,'points')
        xtickangle(0)
    end

end

% Fig 4L
subplot(1,numConds+1,numConds+1)
hold on
for c=1:numConds
    numCells = size(std_by_period{c},4);
    for i=1:numCells
        line([1 2]+(c-1)*2, stdevs_preVsPostDiff{c}(1:2,1,2,i), 'Color',[0 0 0 0.1])
        scatter([1 2]+(c-1)*2, stdevs_preVsPostDiff{c}(1:2,1,2,i), circleSize, 'k','MarkerEdgeAlpha',0.1)
    end
    for o=1:2
        line([o-meanLineWidth o+meanLineWidth]+(c-1)*2, ...
            [mean(stdevs_preVsPostDiff{c}(o,1,2,:),4) mean(stdevs_preVsPostDiff{c}(o,1,2,:),4)], 'Color', [0 0.5 1], 'LineWidth',meanLineThickness)
    end
end
line([0.5 4.5],[0 0],'Color','k')
xlim([0.5 4.5])
xticks(1:4);
xticklabels({'CS+','CS-','CS+','CS-'})
ylabel('\Delta st. dev. (post-pre)')
if forFigure
    set(gca,'units','centimeters','position',[1+7.733,1.5,2.2,2.5])
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    fontsize(gca,6,'points')
end

%% p-value for Fig 4K

% get post minus pre for the learning condition
learning_stdevs_preVsPostDiff = diff(std_by_period{2},1,2);
% (1,:,2,:) vs (2,:,2,:) 
% - 1st index is for CS+ vs CS-
% - 2nd index is a singleton dimension (was pre vs post in diff(std_by_period{2},1,2);
% - 3rd index is the period (2 = during odor)
% - 4th index is the repetitions
pvalue = signrank(squeeze(learning_stdevs_preVsPostDiff(1,:,2,:)), squeeze(learning_stdevs_preVsPostDiff(2,:,2,:)));
fprintf('p-value for Fig 4K: %.5f\n', pvalue); 
