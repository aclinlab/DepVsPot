x = -100:100;
sigmaValues = [2 10 20 50];
numT = length(sigmaValues);

cmap = copper;
figure, hold on
legendText = cell(numT,1);
for ti = 1:numT
    sigma = sigmaValues(ti);
    legendText{ti} = append('\sigma = ',num2str(sigmaValues(ti)));
    color = cmap(ceil(ti/numT*256),:);
    plot(x, normcdf(x/sigma), 'Color', color);
end
line([min(x) max(x)], [0.5 0.5], 'Color', 'k');
line([0 0], [0 1], 'Color', 'k');

xlabel('Avoid - Approach')
ylabel('P(avoid)')
set(gca,'units','centimeters','position',[1,1,2,2])
axis square
set(gca,'TickLength',[.03 .03],'TickDir','out')
fontsize(6,'points')
legend(legendText);