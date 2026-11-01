function [amplitudes, stdevMinValues, rangePeakTimes] = analyzeTraces(input)
%% extract amplitude of first EPSC
startWindow = 1e4;
endWindow = 1.3e4;

% Get the amplitude of each individiual EPSC (ignore previous)
[amplitudes, peakTimes] = min(input(startWindow:endWindow,:),[],1);
rangePeakTimes = range(peakTimes);

% Get the stdev of minimum values within the range during the baseline
% period
baselineStart = 0.5e4;
baselineEnd = 0.8e4;
stdevMinValues = zeros(baselineStart-(baselineEnd-rangePeakTimes)+1,1);
i=1;
for t=baselineStart:(baselineEnd-rangePeakTimes)
    thisWindow = t:(t+rangePeakTimes);
    minValues = min(input(thisWindow,:),[],1);
    stdevMinValues(i) = std(minValues);
    i=i+1;
end
amplitudes = -amplitudes;

end