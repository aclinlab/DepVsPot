load('Fig7G_result.mat');
load('JulianMelissaData.mat')
load('fractions.mat')

numDilutions = size(decisionsTotal,1);
numVoltages = size(decisionsTotal,2);

%% get parameters for deltaXstay habituation
mudeltaXstay_integV = mean(mudeltaXstay_finals,[2 3]);
muA_quad=fitnlm(fractions,mudeltaXstay_integV,'y~b0*x^2+b1',[-1 0.2]);
muA_qcoefs = muA_quad.Coefficients.Estimate;
mudeltaXstay_finals_mean = mean(mudeltaXstay_finals,3);
sum_deltaXstays = sum(mean(mudeltaXstay_finals,3),1); % actual sum
sum_b0_only = sum(fractions.^2*muA_qcoefs(1)); %predicted with b0 only
offsets = (sum_deltaXstays - sum_b0_only) / numDilutions; % b1


colors = [0 0 0;
    0.102, 0.373, 0.478;
    0.739, 0.390, 0.571;
    1, 0.647, 0];
xpoints = 0:0.01:0.26;

markerSize = 3;

%% Fig S4B,C
figure,
subplot(1,2,1)
hold on
plot(fractions,mean(mudeltaXstay_finals,[2 3]),'o','Color','k','MarkerSize',markerSize);
plot(xpoints, muA_qcoefs(1)*xpoints.^2 + muA_qcoefs(2),'-','Color','k')
line([0 0.26],[0 0],'Color','k')
xlim([0 0.26])
ylim([min(mudeltaXstay_finals_mean(:)) max(mudeltaXstay_finals_mean(:))])

ylabel('Change in \mu_{stay}')
xlabel('Fraction KCs responding')
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off


subplot(1,2,2)
hold on
for v=1:numVoltages
    plot(fractions, mean(mudeltaXstay_finals(:,v,:),3),'o--','Color',colors(v,:),'MarkerSize',markerSize)
    plot(xpoints, offsets(v) + muA_qcoefs(1)*xpoints.^2,'-','Color',colors(v,:))
end
xlim([0 0.26])
ylim([min(mudeltaXstay_finals_mean(:)) max(mudeltaXstay_finals_mean(:))])
line([0 0.26],[0 0],'Color','k')


ylabel('Change in \mu_{stay}')
xlabel('Fraction KCs responding')
set(gca,'units','centimeters','position',[4,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off


%% Fig 7G
figure
mergeEnterStay = cat(4,mudeltaXenter_finals,mudeltaXstay_finals);
mergeAll = cat(4,mudeltaXenter_finals,mudeltaXstay_finals,sigmadeltaXenter_finals,sigmadeltaXstay_finals);
mergeSigmas = cat(4,sigmadeltaXenter_finals,sigmadeltaXstay_finals);
commonYlimMu = [min(mergeEnterStay,[],'all'), max(mergeEnterStay,[],'all')];
commonYlimSigma = [0, max(mean(mergeSigmas,3)+std(mergeSigmas,0,3),[],'all')];
commonYlims = [commonYlimMu; commonYlimMu; commonYlimSigma; commonYlimSigma];

faceColors = repmat(colors,1,1,4);
faceColors(:,:,[2 4]) = 1; % stay parameters will be empty circles

ylabels = {'\mu_{enter}','\mu_{stay}','\sigma_{enter}','\sigma_{stay}'};
for i=1:4
    subplot(1,4,i), hold on
    for v=1:numVoltages
        errorbar(fractions, mean(mergeAll(:,v,:,i),3), std(mergeAll(:,v,:,i),0,3), ...
            'o-','Color',colors(v,:),'MarkerFaceColor',faceColors(v,:,i),'MarkerSize',3,'CapSize',3);
    end
    ylabel(['Change in ',ylabels{i}])
    xlabel('Fraction KCs responding')
    % h=legend({'0 V', '15 V', '30 V', '60 V'},'Location','southwest');
    % h.AutoUpdate = 'off';
    line([0 0.26],[0 0],'Color','k')
    xlim([0 0.26])
    ylim(commonYlims(i,:))
    set(gca,'units','centimeters','position',[1+(i-1)*3,1,2,2])
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    axis square
    fontsize(6,'points')
    box off
    
end
