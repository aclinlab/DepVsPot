x = 0:0.1:6;
k=1;
% 1 approach, 2 avoid
% 4 approach, 5 avoid
xvals = [1 2 4 5];

% for Fig 1B, draw the error bars and shading in Illustrator

% for Fig 1C
figure
%% f(x) = x
% subplot(1,3,1);
plotFigure(x, x, xvals)
ylim([min(x) max(x)])

%% log x
figure
% subplot(1,3,2);
plotFigure(x, log(x+k), xvals)


function [] = plotFigure(x, y, xvals)
hold on
plot(x, y, 'k','LineWidth',1);
[~, xis, ~] = intersect(x, xvals);
for i=1:length(xvals)
    if mod(i,2)
        color = [0 0.4470 0.7410]; %blue for approach
    else
        color = [0.9290 0.6940 0.1250]; %orange for avoid
    end
    plot(xvals(i), y(xis(i)), 'o', MarkerEdgeColor=color, MarkerFaceColor=color, MarkerSize=3);
    % vertical line
    line([xvals(i) xvals(i)], [0 y(xis(i))], 'Color',color);
    % horizontal line
    line([0 xvals(i)], [ y(xis(i)) y(xis(i)) ], 'Color',color);
end
xlabel('Input');
ylabel('Output')

set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')

axis square
fontsize(6,'points')
box off
xlim([min(x) max(x)])

end

