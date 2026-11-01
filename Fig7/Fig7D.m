load('fractions.mat');
load('JulianMelissaData.mat');
load('integV.mat');

numDilutions = size(decisionsTotal,1);
numVoltages = size(decisionsTotal,2);


%% naive attraction increases with odor concentration

meanNaiveApp = zeros(numDilutions,1);
semNaiveAv = zeros(numDilutions,1);
ci95NaiveAv = zeros(numDilutions,1);
for d=1:numDilutions
    meanNaiveApp(d) = mean(trainedOdorDecsIntegV{d}(:,2)); % Avoid = 1 - Approach
    semNaiveAv(d) = std(trainedOdorDecsIntegV{d}(:,2))/sqrt(size(trainedOdorDecsIntegV{d},1));
    ci95NaiveAv = ConfInt95(trainedOdorDecsIntegV{d}(:,2),1);
end
figure
hold on
errorbar(fractions,meanNaiveApp,ci95NaiveAv,'o-','Color','k','MarkerSize',3,'CapSize',3,'MarkerFaceColor','k')
plot(fractions,meanNaiveApp,'Color','k')
xl = [0 0.26];
xlim(xl);
line(xl,[0.5 0.5],'Color','k');
xlabel('Fraction KCs responding');
ylabel('p_{obs}(approach)')
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out','Color','none')
axis square
fontsize(6,'points')
fontname('Arial')
box off

%% P(enter), P(stay) naive
for d=1:numDilutions
    meanPenter_target_integV(d,:) = mean(decsTrainedOdorEnterIntegV{d},1,'omitnan');
    n_enter = sum(~isnan(decsTrainedOdorEnterIntegV{d}),1);
    semPenter_target_integV(d,:) = std(decsTrainedOdorEnterIntegV{d},0,1,'omitnan')./sqrt(n_enter);
    ci95Penter_target_integV(d,:) = semPenter_target_integV(d,:).*tinv(0.975, n_enter-1);
    meanPstay_target_integV(d,:) = mean(decsTrainedOdorStayIntegV{d},1,'omitnan');
    n_stay = sum(~isnan(decsTrainedOdorStayIntegV{d}),1);
    semPstay_target_integV(d,:) = std(decsTrainedOdorStayIntegV{d},0,1,'omitnan')./sqrt(n_stay);
    ci95Pstay_target_integV(d,:) = semPstay_target_integV(d,:).*tinv(0.975, n_stay-1);
end

figure,hold on
if reverseLeave
    errorbar(fractions,1-meanPenter_target_integV(:,2),ci95Penter_target_integV(:,2),'o-','MarkerSize',3,'Color','k','MarkerFaceColor','k','CapSize',3);
    errorbar(fractions,1-meanPstay_target_integV(:,2),ci95Pstay_target_integV(:,2),'o-','MarkerSize',3,'Color','k','MarkerEdgeColor','k','MarkerFaceColor','w','CapSize',3);
else
    errorbar(fractions,meanPenter_target_integV(:,2),ci95Penter_target_integV(:,2),'o-','MarkerSize',3,'Color','k','MarkerFaceColor','k','CapSize',3);
    errorbar(fractions,meanPstay_target_integV(:,2),ci95Pstay_target_integV(:,2),'o-','MarkerSize',3,'Color','k','MarkerEdgeColor','k','MarkerFaceColor','w','CapSize',3);
end
xl = [0 0.26];
xlim(xl);
line(xl,[0.5 0.5],'Color','k');
xlabel('Fraction KCs responding');
if reverseLeave
    ylabel('p_{obs}(reverse) or p_{obs}(leave)')
else
    ylabel('p_{obs}(enter) or p_{obs}(stay)')
end
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out','Color','none')
axis square
fontsize(6,'points')
fontname('Arial')

box off