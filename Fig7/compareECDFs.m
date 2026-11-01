function [] = compareECDFs(realData, simData)
hold on
bins = 0.005:0.01:0.995;
realECDF = interpECDF(realData,bins);
stairs([bins-0.005, 1],[realECDF,1],'r');
numTrials = size(simData,2);
simHistValues = zeros(numTrials,length(bins));
for t=1:numTrials
    try
        [simHistValues(t,:),~] = interpECDF(simData(:,t),bins);
    catch
        error('here we are')
    end
end
errorRange = prctile(simHistValues,[2.5 97.5],1);
meanSimHistValues = mean(simHistValues,1,'omitnan');
stairs([bins-0.005, 1],[meanSimHistValues,1],'k');
binEdges = [0 repelem((0.01:0.01:0.99),2) 1];
fillx = [binEdges, fliplr(binEdges)];
filly = repelem([errorRange(1,:) fliplr(errorRange(2,:))],2);
fill(fillx, filly,'k','FaceAlpha',0.2,'Edgecolor','none')
ylim([0 1])

end