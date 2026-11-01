S=load('MBmodels.mat');
n=2000;
odors = 176;
numTrials=30;
m=24;
% use the activity-independent adaptation model - first model instance
KCresponse = reshape(calculateKCresponse(S.MBmodels_all{2,1},reshape(S.allTestedOdors{1},m,[]),false),n,odors,numTrials);
KCresponse = KCresponse(:,:,1); % sum(KCresponse(:,:,1:15),3); % sum over training trials
KCresponse = KCresponse/mean(KCresponse(KCresponse>0),'all');
Ta = 0.0215*n*mean(KCresponse(KCresponse>0),'all');
Tm = 0.1710;
mu_enter_0 = 1.4436 - 0.2885;

% %% find the exponential distribution fit for each odor
% 
% params = nan(odors,2);
% figure, hold on
% for i=1:odors
%     z = KCresponse(:,i);
%     z = z(z>0);
%     if length(z)>10
%         z = z/mean(z);
%         zs = sort(z);  F = (1:numel(zs))'/numel(zs);    % empirical cumulative distribution
%         lowEnd = zs>0;% < 0.3;                              % focus on the weak responses
%         loglog(zs(lowEnd), F(lowEnd));  xlabel('response x');  ylabel('fraction below x')
%         p = polyfit(log10(zs(lowEnd)), log10(F(lowEnd)), 1);
%         params(i,1) = p(1); %(1 = flat near zero, as for an exponential)
%         params(i,2) = 10^p(2); %p0
%     end
%     % fprintf('k = %.2f (1 = flat near zero, as for an exponential)\n', p(1))
%     % fprintf('p0 = %.2f\n', 10^p(2))                 % if k is about 1: density at zero
% end

%% per pair find eta-star

act   = KCresponse > 0;
nOd   = size(KCresponse, 2);
pairs = nchoosek(1:nOd, 2);
nP    = size(pairs, 1);
etaStar = zeros(nP,1);
logLRs = -2:0.2:2;
nLRs = length(logLRs);
learningRates = 10.^logLRs;
learning = zeros(nP,nLRs,2); % num pairs x num learning rates x dep/pot
DeltaXs = zeros(nP,nLRs,2);
% loop through all possible pairs

parfor p = 1:nP
    % p
    idx1 = pairs(p,1);
    idx2 = pairs(p,2);
    idx_overlap = act(:, idx1) & act(:, idx2); % indices of overlapping KCs
    idx_odor1_only = act(:, idx1) & ~act(:, idx2); % KCs that respond only to odor 1
    etaStar(p) = etaStarNewton(KCresponse(idx_odor1_only, idx1,1), ...
        KCresponse(idx_overlap, idx1,1), ...
        KCresponse(idx_overlap, idx2,1) );

    thisPair = zeros(nLRs,2);
    thisPairX = zeros(nLRs,2);
    for i=1:nLRs

        learningRate = learningRates(i);

        % % column 1 is Avoid, column 2 is Approach
        wDep = exp(-learningRate*KCresponse(:,[idx2 idx1]));% numKCs x Av/App - depress wrong actions by exponential decay. column 2 of KC activity is rewarded so depress avoid (column 1)
        wPot = 1-exp(-learningRate*KCresponse(:,[idx1 idx2]));

        DepO = wDep' * KCresponse(:,[idx1 idx2]); % stands for Depression Output
        PotO = wPot' * KCresponse(:,[idx1 idx2]);

        DeltaXavoidDep = (DepO(1) - DepO(2)) / sqrt( (Ta+Tm*DepO(1))^2 + (Ta+Tm*DepO(2))^2 );
        DeltaXavoidPot = (PotO(1) - PotO(2)) / sqrt( (Ta+Tm*PotO(1))^2 + (Ta+Tm*PotO(2))^2 )
        PavoidDep = 1-normcdf(mu_enter_0 - DeltaXavoidDep );
        PavoidPot = 1-normcdf(mu_enter_0 - DeltaXavoidPot );

        thisPair(i,:) = [PavoidDep; PavoidPot];
        thisPairX(i,:) = [DeltaXavoidDep; DeltaXavoidPot];
    end
    learning(p,:,:) = thisPair;
    DeltaXs(p,:,:) = thisPairX;

end

%% get learning at etaStar as fraction of final learning

learningAtEtaStar = zeros(nP,1);
learningAtEtaStarFrac = zeros(nP,1);
DeltaXsAtEtaStar = zeros(nP,1);
DeltaXsAtEtaStarFrac = zeros(nP,1);
DeltaWAtEtaStar = zeros(nP,1);
lowerBoundEtaStar = zeros(nP,1);
upperBoundEtaStarWithAssumptions = zeros(nP,1);
spearmanCorrOverlap = zeros(nP,1);
pearsonCorrOverlap = zeros(nP,1);
nOverlap = zeros(nP,1);
odor1odor2ratio = zeros(nP,1);
overlapodor1ratio = zeros(nP,1);
Rs = zeros(nP,1);

c=1;
for p=1:nP
    idx1 = pairs(p,1);
    idx2 = pairs(p,2);
    idx_overlap = act(:, idx1) & act(:, idx2); % indices of overlapping KCs
    nOverlap(p) = nnz(idx_overlap);
    idx_odor1_only = act(:, idx1) & ~act(:, idx2); % KCs that respond only to odor 1
    S = sum(KCresponse(idx_odor1_only,idx1)) + sum(KCresponse(idx_overlap,idx1));
    Q = sum(KCresponse(idx_odor1_only,idx1).^2) + sum(KCresponse(idx_overlap,idx1).*(KCresponse(idx_overlap,idx1) + KCresponse(idx_overlap,idx2)));
    lowerBoundEtaStar(p) = S/Q;

    R = (nnz(idx_overlap)+nnz(idx_odor1_only))/sum(KCresponse(idx_overlap,idx1));
    Rs(p) = R;
    basedOnExpDist = (1+sqrt(1+4*R))/2;
    upperBoundEtaStarWithAssumptions(p) = min(R, basedOnExpDist); %(nnz(idx_overlap)+nnz(idx_odor1_only))/sum(KCresponse(idx_overlap,idx1));
    % upperBoundEtaStarWithAssumptions(p) = (nnz(idx_overlap)+nnz(idx_odor1_only))/nnz(idx_overlap);

    odor1odor2ratio(p) = mean(1-exp(-basedOnExpDist*KCresponse(idx_overlap,idx2))) / mean(1-exp(-basedOnExpDist*KCresponse(idx_overlap|idx_odor1_only,idx1)));
    overlapodor1ratio(p) = mean(KCresponse(idx_overlap,idx1)) / mean(KCresponse(idx_odor1_only,idx1));

    if nnz(idx_overlap)
        spearmanCorrOverlap(p) = corr(KCresponse(idx_overlap,idx1), KCresponse(idx_overlap,idx2), 'Type', 'Spearman');
        pearsonCorrOverlap(p) = corr(KCresponse(idx_overlap,idx1), 1-exp(-basedOnExpDist*KCresponse(idx_overlap,idx2)), 'Type', 'Pearson');
    else
        spearmanCorrOverlap(p) = nan;
        pearsonCorrOverlap(p) = nan;
    end
    [~,closestIndex] = min(abs(etaStar(p)-learningRates));
    learningAtEtaStar(p) = learning(p,closestIndex,1);
    learningAtEtaStarFrac(p) = learningAtEtaStar(p)/max(learning(p,end,1));
    DeltaXsAtEtaStar(p) = DeltaXs(p,closestIndex,1);
    DeltaXsAtEtaStarFrac(p) = DeltaXsAtEtaStar(p)/DeltaXs(p,end,1);
    DeltaW = 1-exp(-etaStar(p)*KCresponse(:,pairs(p,1)));
    DeltaWAtEtaStar(p) = mean(DeltaW(DeltaW>0));
end

% not in the manuscript- but it shows which assumptions are violated for
% the points above the estimated upper bound
rgbcolors = zeros(nP,3);
rgbcolors(:,1) = spearmanCorrOverlap>0;
rgbcolors(:,2) = overlapodor1ratio >= 1;
rgbcolors(:,3) = odor1odor2ratio >= 1;

figure,scatter(log10(upperBoundEtaStarWithAssumptions),log10(etaStar),5,rgbcolors*255); %,spearmanCorrOverlap>=0 & odor1odor2ratio>=1 & overlapodor1ratio>=1)
hold on, line([-1 3], [-1 3])
set(gca,'Color','k')


%% get overlap
% correlations = zeros(nP,1);
% angDists = zeros(nP,1);
% fracOverlapOdor1 = zeros(nP,1);
fracOverlap = zeros(nP,1);
for p=1:nP
    idx1 = pairs(p,1);
    idx2 = pairs(p,2);
    % correlations(p) = corr(KCresponse(:,idx1), KCresponse(:,idx2));
    % angDists(p) = pdist2(KCresponse(:,idx1)', KCresponse(:,idx2)','cosine');
    % fracOverlapOdor1(p) = nnz(act(:, idx1) & act(:, idx2)) / nnz(act(:, idx1));
    fracOverlap(p) = nnz(act(:, idx1) & act(:, idx2)) / n;
end



%% make the figure

figure('Position',[100 100 600 300])

% cdf of eta-star
ax(1) = subplot(2,4,5);
zs = sort(log10(etaStar));  F = (1:numel(zs))'/numel(zs);
plot(zs,F,'Color','k');
xlabel('log_{10}(\eta*)')
ylabel('Cumulative distribution')
xlim([zs(1) zs(end)])
ylim([0 1]);
line(0.8*[1 1],[0 1],'Color','k')
box off

% pdf of non-zero KC responses compared to exponential
ax(2) = subplot(2,4,6);
hold on
x  = KCresponse(KCresponse>0); 
xs = linspace(0, max(x), 300);
pd = fitdist(x, 'Exponential');
histogram(x, 'Normalization', 'pdf','FaceColor','k','EdgeColor','none'); 
plot(xs, pdf(pd, xs), 'r', 'LineWidth', 0.5);
xlabel('KC response amplitude'); 
ylabel('Density');
box off

% R vs etaStar
ax(3) = subplot(2,4,2);
hold on
passAssumptions = (pearsonCorrOverlap>=0)&(odor1odor2ratio>=1);
scatter(log10(Rs(~passAssumptions)),log10(etaStar(~passAssumptions)),3,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1)
scatter(log10(Rs(passAssumptions)),log10(etaStar(passAssumptions)),3,'MarkerEdgeColor','r','MarkerEdgeAlpha',0.1)
hold on, line([-0.6 6], [-0.6 6],'Color','c')
x=10.^(-1:0.1:6); %[10^-1, 10^6];
p0=0.9;
line(log10(x),log10((p0+sqrt(p0^2+4*p0*x))/2),'Color','g')
ylim([-0.6 6])
xlim([-0.6 6])
line([-0.6 6],0.8*[1 1],'Color','k')
xlabel('log_{10}(R)');
ylabel('log_{10}(\eta*)')

% R vs etaStar colored by fraction overlap
ax(4) = subplot(2,4,3);
hold on
scatter(log10(Rs),log10(etaStar),3,fracOverlap,'MarkerEdgeAlpha',0.1)
hold on, line([-0.6 6], [-0.6 6],'Color','c')
x=10.^(-1:0.1:6); %[10^-1, 10^6];
p0=0.9;
line(log10(x),log10((p0+sqrt(p0^2+4*p0*x))/2),'Color','g')
ylim([-0.6 6])
xlim([-0.6 6])
line([-0.6 6],0.8*[1 1],'Color','k')
xlabel('log_{10}(R)');
ylabel('log_{10}(\eta*)')

meanPerOdorDivByTotalMean = sum(KCresponse,1)./sum(KCresponse>0,1)/mean(KCresponse(KCresponse>0));
ax(5) = subplot(2,4,7);
scatter(sum(KCresponse>0,1)/n,meanPerOdorDivByTotalMean,3)
xl=xlim;
yl=ylim;
line(xl,[1 1],'Color','k')
line(0.005*[1 1],yl,'Color','k')
xlabel('Fraction KCs activated')
ylabel({'Mean response ampl. per odor','Mean response ampl. across odors'})

% avg learning for different LRs
ax(6) = subplot(2,4,4);
plot(logLRs,squeeze(mean(learning,1)))
line(0.8*[1 1],[0.1 0.6],'Color','k')
line(logLRs([1 end]),normcdf(-mu_enter_0)*[1 1],'Color','k')
ylim([0.1 0.6])
box off
xlabel('log_{10}(\eta)')
ylabel('P(reverse), Odor A')

% Dep-Pot at "complete" learning
ax(7) = subplot(2,4,8);
scatter(fracOverlap,-diff(squeeze(learning(:,15,:)),1,2),3,sum(KCresponse(:,pairs(:,1))>0,1)/n)
xlabel('Overlap')
ylabel({'P(reverse), Odor A','Difference (Dep - Pot)'})


arrayfun(@(a) set(a, 'TickDir', 'out'), ax);
arrayfun(@(a) axis(a, 'square'), ax);
for i=[3 4 7]
    rasterizeScatterLayer(ax(i),600);
end
for i=1:length(ax)
    set(ax(i),'Color','none')
    box(ax(i),'off')
    set(ax(i),'units','centimeters','FontSize',6,'TickLength',[.03 .03])
    ax(i).Position(3:4) = [3 3];
end

