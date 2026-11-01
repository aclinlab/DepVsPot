%% Data explanation from Toshihide Hige:

% Each matrix corresponds to data from a single cell.
% n-th column corresponds to n-th trial.
% The columns 5*(n-1)+1 to 5*n correspond to 5 trials from the n-th minute. So first 15 columns are pre-pairing trials.
% The first trials of the minute (i.e. 5*(n-1)+1) are single-pulse trials, and the rest are paired-pulse trials. So for your purpose, you can use the first pulse of any trials.
% Sampling rate is 10 kHz.
% The data are leak-subtracted at the baseline period and low-pass filtered.
% 

%% load data

load('YamadaFig2C.mat') % DA paired
load('YamadaFig5D.mat') % BAY paired
BAY_Pairing_cell1Fig5D = BAY_Pairing_cell1;
BAY_Pairing_cell2Fig5D = BAY_Pairing_cell2;
BAY_Pairing_cell3Fig5D = BAY_Pairing_cell3;
BAY_Pairing_cell4Fig5D = BAY_Pairing_cell4;
BAY_Pairing_cell5Fig5D = BAY_Pairing_cell5;
load('YamadaFig2K.mat') % DA paired + DMSO
load('YamadaFig3F.mat') % low forskolin paired
load('YamadaFig8D.mat') % a/b BAY paired

load('YamadaFig2E.mat') % DA only
load('YamadaFig2G.mat') % LED only
load('YamadaFig3B.mat') % forskolin only
load('YamadaFig3H.mat') % DMSO paired

numConds = 2;
allData = cell(numConds,1); % 1. control, 2. depression

% control
allData{1} = cat(3, DA_Only_cell1,...
  DA_Only_cell2,...
  DA_Only_cell3,...
  DMSO_pairing_cell1,... % lower case pairing is Fig 3H DMSO as neg ctrl for forskolin
  DMSO_pairing_cell2,...
  DMSO_pairing_cell3,...
  Fsk_20uM_Only_cell1,...
  Fsk_20uM_Only_cell2,...
  Fsk_20uM_Only_cell3,...
  Fsk_20uM_Only_cell4,...
  Fsk_20uM_Only_cell5,...
  Fsk_20uM_Only_cell6,...
  LED_Only_cell1,...
  LED_Only_cell2,...
  LED_Only_cell3);

% depression
allData{2} = cat(3,DA_Pairing_cell1,...
    DA_Pairing_cell2,...
    DA_Pairing_cell3,...
    DA_Pairing_cell4,...
    DA_Pairing_cell5,...
    DA_Pairing_cell6,...
    DMSO_Pairing_cell1,... % upper case Pairing is Fig 2K DMSO as neg ctrl for SCH
    DMSO_Pairing_cell2,...
    DMSO_Pairing_cell3,...
    DMSO_Pairing_cell4,...
    DMSO_Pairing_cell5,...
    DMSO_Pairing_cell6,...
    Fsk_20uM_Pairing_cell1,...
    Fsk_20uM_Pairing_cell2,...
    Fsk_20uM_Pairing_cell3,...
    Fsk_20uM_Pairing_cell4,...
    Fsk_20uM_Pairing_cell5,...
    Fsk_20uM_Pairing_cell6);

%% analyze data

[stdEPSCAmp, stdBaseline, meanEPSCAmp, rangeForStdev, allPreAmplitudes, allPostAmplitudes, allSmoothPreAmps, allSmoothPostAmps] = deal(cell(numConds,1));

for i=1:length(allData)
    numDataSets = size(allData{i},3);
    [stdEPSCAmp{i}, stdBaseline{i}, meanEPSCAmp{i}, rangeForStdev{i}] = deal(zeros(numDataSets,2));
    [allPreAmplitudes{i}, allPostAmplitudes{i}, allSmoothPreAmps{i}, allSmoothPostAmps{i}] = deal(cell(numDataSets,1));

    for d=1:numDataSets %d for data
        d
        data = allData{i}(:,:,d);
        totalTime = size(data,1);

        % pre-pairing
        prepairing = data(:,1:15);
        [amplitudes, stdevMinValues, rangePeakTimes] = analyzeTraces(prepairing);
        smoothAmps = malowess(1:length(amplitudes),amplitudes,Span=min(length(amplitudes),25));
        deviations = amplitudes - smoothAmps;
        stdEPSCAmp{i}(d,1) = sqrt(sum(deviations.^2)/(length(deviations)-1)); % stdev not from the total mean but from the smoothed trend std(amplitudes); %
        meanEPSCAmp{i}(d,1) = mean(amplitudes);
        stdBaseline{i}(d,1) = mean(stdevMinValues);
        rangeForStdev{i}(d,1) = rangePeakTimes;
        allPreAmplitudes{i}{d} = amplitudes;
        allSmoothPreAmps{i}{d} = smoothAmps;

        % post-pairing
        postpairing = data(:,21:end);
        [amplitudes, stdevMinValues, rangePeakTimes] = analyzeTraces(postpairing);
        smoothAmps = malowess(1:length(amplitudes),amplitudes,Span=min(length(amplitudes),25));
        deviations = amplitudes - smoothAmps;
        stdEPSCAmp{i}(d,2) = sqrt(sum(deviations.^2)/(length(deviations)-1)); 
        meanEPSCAmp{i}(d,2) = mean(amplitudes);
        stdBaseline{i}(d,2) = mean(stdevMinValues);
        rangeForStdev{i}(d,2) = rangePeakTimes;
        allPostAmplitudes{i}{d} = amplitudes;
        allSmoothPostAmps{i}{d} = smoothAmps;

    end
end

%% Supp Fig: plot EPSC amplitudes by trial
for i=1:numConds
    figure
    figHandle = gcf;
    set(figHandle,'units','centimeters')
    figHandle.Position(3:4) = [21 9];
    numDataSets = size(allData{i},3);
    for d=1:numDataSets
        subplot(3,6,d)
        hold on
        scatter(1:15,allPreAmplitudes{i}{d},3,'MarkerEdgeAlpha',0.5)
        scatter(21:90,allPostAmplitudes{i}{d},3,'MarkerEdgeAlpha',0.5)
        plot(1:15,allSmoothPreAmps{i}{d},'k')
        plot(21:90,allSmoothPostAmps{i}{d},'k')
        xlim([0 90])
        xticks(0:30:90);

        yl = ylim;
        ylim([0 yl(2)])
        xlabel('Trial')
        ylabel('EPSC amplitude (pA)')
        axHandle = gca;
        set(axHandle,'units','centimeters')
        axHandle.Position(3:4) = [1.8 1.8];
        set(gca,'TickLength',[.03 .03],'TickDir','out')
        axis square
        fontsize(6,'points')
        box off
    end
end


%% Fig 3B: example traces

figure,hold on
plot(allData{2}(100:100:end,1:15,16),'k') % downsample traces to allow exporting vector graphics
plot(allData{2}(100:100:end,16:end,16),'b') % downsample traces to allow exporting vector graphics
xlim([0.5e2 1.3e2])
% The sample rate is 10 kHz
% by down-sampling by a factor of 100, the x-axis is in units of 10 ms. So
% if the tick marks are 20, that's 200 ms
xlabel('Time (10xms)')
ylabel('current (pA)')
axHandle = gca;
set(axHandle,'units','centimeters')
axHandle.Position(3:4) = [2.2 2.2];
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off

%% Fig 3C: plot EPSC amplitudes by trial

whichToPlot = [10 16];
figure
figHandle = gcf;
set(figHandle,'units','centimeters')
figHandle.Position(3:4) = [9 6];
for i=1:2 % control and depression
    numDataSets = size(allData{i},3);
    d = whichToPlot(i);
    subplot(1,2,i)
    hold on
    scatter(1:15,allPreAmplitudes{i}{d},3,'k','MarkerEdgeAlpha',0.5)
    scatter(21:90,allPostAmplitudes{i}{d},3,'b','MarkerEdgeAlpha',0.5)
    plot(1:15,allSmoothPreAmps{i}{d},'k')
    plot(21:90,allSmoothPostAmps{i}{d},'b')
    xlim([0 90])
    xticks(0:30:90);

    yl = ylim;
    ylim([0 yl(2)])
    xlabel('Trial')
    ylabel('EPSC amplitude (pA)')
    axHandle = gca;
    set(axHandle,'units','centimeters')
    axHandle.Position(3:4) = [2.5 2.5];
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(6,'points')
    box off

end



%% Fig 3D,E lines
depStart = 1;
depEnd = 18;
potStart = 19;
potEnd = 23;
condColors = {'black', 'red', 'green'};
for i=1:numConds
    normVars{i} = stdEPSCAmp{i}.^2-stdBaseline{i}.^2;
    diffNormVars{i} = normVars{i}(:,2) - normVars{i}(:,1);
    normStd{i} = sqrt(normVars{i});
    diffNormStd{i} = normStd{i}(:,2)-normStd{i}(:,1);
    VarOverMean{i} = normVars{i}./meanEPSCAmp{i};
    figure
    hold on;
    % before
    scatter(meanEPSCAmp{i}(:,1),normStd{i}(:,1),8, condColors{i});
    % after
    scatter(meanEPSCAmp{i}(:,2),normStd{i}(:,2),8, condColors{i}, 'filled');
    % connecting lines
    numDataSets = size(allData{i},3);
    for d=1:numDataSets
        line([meanEPSCAmp{i}(d,:)], normStd{i}(d,:),'Color',condColors{i});
    end
    xlabel('Mean EPSC amplitude (pA)')
    ylabel({'Baseline-corrected', 'standard deviation (pA)'})
    axHandle = gca;
    set(axHandle,'units','centimeters')
    axHandle.Position(3:4) = [3.33 3.33];
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(6,'points')
    box off

end

%% diffNormStd gives the data plotted in Fig 3F
disp('Mann-Whitney test p-value for ∆ baseline-corrected std, ctrl vs depression')
ranksum(diffNormStd{1},diffNormStd{2},'method', 'exact');

