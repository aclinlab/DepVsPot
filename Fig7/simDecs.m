function result = simDecs(numDecsAllFlies, PenterAllFlies, PstayAllFlies)
% result = simDecs(numDecsAllFlies, PenterAllFlies, PstayAllFlies)
% Simulates decision fractions for N flies
% inputs are all N x 1 arrays where N = numFlies
%  numDecsAllFlies is the number of decisions of each fly made
%  PenterAllFlies is P(enter) [probability of entering the CS+ side when
%   starting on the other side]. Note this is Nx1, every fly might have a
%   different P(enter)
%  PstayAllFlies - same as PenterAllFlies, but for P(stay) [probability of
%   staying on the CS+ side when starting on that side]
% Note: if the input arrays are not the same size:
%   if numDecsAllFlies is smaller, the function will not give an error but
%   it will ignore the flies in PenterAllFlies and PstayAllFlies that are
%   at indices "after" the end of numDecsAllFlies
%   If numDecsAllFlies is bigger, it will crash with an "out of range" error
% The function simulates a single random trial, where at each
% decision, each fly "decides" what to do based on the specified P(enter)
% and P(stay). 
% result will be N x 3 - 
% each row is one fly's fractions of decisions to approach, stay, enter
% column 1: approach 
% column 2: stay 
% column 3: enter

numFlies = length(numDecsAllFlies);
result = zeros(numFlies,3);
% loop over all flies
for i=1:numFlies
    numDecs = numDecsAllFlies(i); % get the number of decisions
    % onTrainedOdor keeps track of whether the fly is on CS+ side or not: 1 for trained odor, 0 for other side
    onTrainedOdor = logical(round(rand(1))); % start on a random side
    % initialise counts
    decsFor = 0;
    decsStay = 0;
    decsEnter = 0;
    decsStartIn = 0;
    decsStartOut = 0;
    
    % get P(enter) and P(stay)
    Penter = PenterAllFlies(i);
    Pstay = PstayAllFlies(i);

    % simulate the trial where the fly makes a certain number of decisions
    for j = 1:numDecs
        if onTrainedOdor % on trained odor side
            stay = rand(1) < Pstay; % stay will be 1 with probability P(stay)
            onTrainedOdor = stay; % stay on trained odor side
            % update the counters
            decsFor = decsFor + stay;
            decsStay = decsStay + stay;
            decsStartIn = decsStartIn + 1;
        else % on other side
            % This naturally simulates a truncated
            % Gaussian distr of PEnter bc rand(1) is on
            % (0,1)
            enter = rand(1) < Penter; % enter will be 1 with probability P(enter)
            onTrainedOdor = enter; % go from other side to trained side
            % update the counters
            decsFor = decsFor + enter;
            decsEnter = decsEnter + enter;
            decsStartOut = decsStartOut + 1;
        end
    end
    % calculate the fractions
    result(i,:) = [decsFor/numDecs, decsStay/decsStartIn, decsEnter/decsStartOut];
end
end