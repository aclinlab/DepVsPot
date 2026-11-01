colors = [0 0 0;
    0.102, 0.373, 0.478;
    0.739, 0.390, 0.571;
    1, 0.647, 0];


load('fractions.mat'); % fractions of KCs activated by each concentration

load('JulianMelissaData.mat'); % behaviour data - see Readme

numDilutions = size(decisionsTotal,1);
numVoltages = size(decisionsTotal,2);
%% decision scores

[mean_deltaPapp, sem_deltaPapp] = deal(zeros(numDilutions,numVoltages));
for d=1:numDilutions
    for v=1:numVoltages
        deltaPav = trainedOdorDecs{d,v}(:,6) - trainedOdorDecs{d,v}(:,2);
        mean_deltaPapp(d,v) = mean(deltaPav);
        sem_deltaPapp(d,v) = std(deltaPav)/sqrt(length(deltaPav));
    end
end
figure, hold on

for v=1:numVoltages
    errorbar(fractions,mean_deltaPapp(:,v),sem_deltaPapp(:,v),'o-','Color',...
        colors(v,:),'MarkerSize',3,'CapSize',3,'MarkerFaceColor',colors(v,:),'MarkerEdgeColor','none');

end
line([0 0.26],[0 0],'Color','k')
xlim([0 0.26])
set(gca,'units','centimeters','position',[1,1,2.5,2.5])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off
xlabel('Fraction KCs responding')
ylabel('\Deltap_{obs}(approach CS+)')
