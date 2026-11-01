aValues = 1:100;
bValues = 1:100;
numA = length(aValues);
numB = length(bValues);
numTrials = 100;
figure
figHandle = gcf;
set(figHandle,'units','centimeters')
figHandle.Position(3:4) = [14 11];

numCols = 3;
%% additive noise
T = 20*sqrt(2); %stdev of each output is 20
result = zeros(numA, numB);
for ai = 1:numA
    a = aValues(ai);
    for bi = 1:numB
        b = bValues(bi);
        result(ai,bi) = normcdf((a-b)/T);
    end
end
subplot(1,numCols,1);
plotHeatmap(result);

%% Ta, Tm
Ta = 1.2;
Tm = 0.5;
result = zeros(numA, numB);
for ai = 1:numA
    a = aValues(ai);
    for bi = 1:numB
        b = bValues(bi);
        result(ai,bi) = normcdf( (a-b)/(sqrt( (Ta+Tm*a)^2 + (Ta+Tm*b)^2 ) ) ); 
    end
end
subplot(1,numCols,2);
plotHeatmap(result);

%% log transfer function
T = .625*sqrt(2); %stdev of each output is 0.625
result = zeros(numA, numB);
for ai = 1:numA
    a = aValues(ai);
    for bi = 1:numB
        b = bValues(bi);
        result(ai,bi) = normcdf((log(a+1)-log(b+1))/T );
    end
end
subplot(1,numCols,3);
plotHeatmap(result);


%% functions
function [] = plotHeatmap(result)
imagesc(flipud(result'));
hold on
contour(smooth2a(flipud(result'),3,3),0.1:0.1:0.9,'Color','k');
axis equal, axis tight
axHandle = gca;
set(axHandle,'XTick',[]);
set(axHandle,'YTick',[]);
set(axHandle,'units','centimeters'),
axHandle.Position(3:4) = [2,2];
% set(axHandle,'units','normalized'),
end
