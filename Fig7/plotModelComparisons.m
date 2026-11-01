function [] = plotModelComparisons(data,model,modelError,figTitle,varargin)

if nargin>4
    forFigure = varargin{1};
else
    forFigure = false;
end

load('fractions.mat','fractions');


colors = [0 0 0;
    0.102, 0.373, 0.478;
    0.739, 0.390, 0.571;
    1, 0.647, 0];


numTmConds = size(model,3);
numVoltages = size(data,2);
figure
if ~forFigure
    sgtitle(figTitle)
end
min_ylim = Inf;
max_ylim = -Inf;
if forFigure
    markerSize = 3;
    capSize = 3;
else
    markerSize = 6;
    capSize = 6;
end
for TmCond=1:numTmConds
    axHandles(TmCond) = subplot(1,numTmConds,TmCond);
    switch TmCond
        case 1
            title('Tm=0')
        case 2
            title('fit Tm')
    end
    hold on
    for v=1:numVoltages
        plot(fractions, data(:,v),'o--','Color',colors(v,:),'MarkerSize',markerSize);
    end
    for v=1:numVoltages
        errorbar(fractions, model(:,v,TmCond), model(:,v,TmCond)-modelError(:,v,1,TmCond),...
            modelError(:,v,2,TmCond)-model(:,v,TmCond),"-x",'Color',colors(v,:),'MarkerSize',markerSize,'CapSize',capSize);
    end

    ylabel(figTitle)
    xlabel('Fraction KCs responding')
    % h=legend({'0 V', '15 V', '30 V', '60 V', '0 V sim', '15 V sim', '30 V sim', '60 V sim', },'Location','southwest');
    xl = xlim;
    xl(2) = 0.26;
    xlim([0 xl(2)])
    yl = ylim;
    if yl(1) < min_ylim
        min_ylim = yl(1);
    end
    if yl(2) > max_ylim
        max_ylim = yl(2);
    end
    
    if forFigure
        disp(forFigure)
        set(axHandles(TmCond),'units','centimeters','position',[1+3*(TmCond-1),1,2,2])
        set(axHandles(TmCond),'TickLength',[.03 .03],'TickDir','out')
        axis square
        fontsize(axHandles(TmCond),6,'points')
    end
end
for TmCond=1:numTmConds
    set(axHandles(TmCond),'YLim', [min_ylim max_ylim]);
end

end