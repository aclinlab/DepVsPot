
S = load('naiveCombined10Reps.mat');
out = S.out;
% 3 x 10. 
% rows is the models, columns is repetitions
% row 1 = fit model to synthetic data with no variability (sigmaZero)
% row 2 = "                                equal sigma across dilutions (slopeZero)
% row 3 = real data

[numConds, numReps] = size(out);
N = 2000;

modelMakeup = zeros(3, numConds, numReps); % models
meanThetas = zeros(4,3,numConds,numReps); % theta, models
hdi95s = zeros(4,2,3,numConds,numReps); % theta, lower/upper, models
for c=1:numConds
    for p=1:numReps
        for m=1:3
            m_indices = (out(c,p).theta(:,5)==m);
            modelMakeup(m,c,p) = sum(m_indices)/N;
            meanThetas(:,m,c,p) = out(c,p).weights(m_indices)' * out(c,p).theta(m_indices,1:4);
            for d=1:4
                hdi95s(d,:,m,c,p) = hdi95_weighted(out(c,p).theta(m_indices,d), out(c,p).weights(m_indices));
            end
        end
    end
    fprintf('Model %d: %f correct\n', c, mean(modelMakeup(c,c,:)))
end
[~,winningModel] = max(modelMakeup,[],1);
winningModel = squeeze(winningModel); % numConds x numPriors


%% Fig 7F

figure
c=3; % which condition - real data

forPub = true;
verticalAlignment = true;
reverseLeave = false;
if reverseLeave
    xlabels = {'\mu_{reverse}','\sigma_{reverse}','m_{\mu,leave}','m_{\sigma,leave}'};
else
    xlabels = {'\mu_{enter}','\sigma_{enter}','m_{\mu,stay}','m_{\sigma,stay}'};
end
thetas_all_linear = [];
weights_all_linear = [];

for p=1:size(out,2)
    indices = out(c,p).theta(:,5)==winningModel(c);
    thetas_all_linear = [thetas_all_linear; out(c,p).theta(indices,:)];
    weights_all_linear = [weights_all_linear; out(c,p).weights(indices)];
end
weights_all_linear = weights_all_linear/sum(weights_all_linear);

for d=1:4
    subplot(1,4,d), hold on
    signflip = 1;
    if reverseLeave
        if d==1 || d==3
            signflip = -1;
        end
    end
    % separate trace for each run
    for p=1:size(out,2)
        indices = out(c,p).theta(:,5)==winningModel(c);
        [f,xi] = ksdensity(signflip*out(c,p).theta(indices,d),'Weights',out(c,p).weights(indices),'BoundaryCorrection','reflection');
        plot(xi,f,'Color',[0 0 0 0.2])
    end
    % all particles over all runs
    [f,xi] = ksdensity(thetas_all_linear(:,d),'Weights',weights_all_linear,'BoundaryCorrection','reflection');
    plot(xi,f,'Color','k')

    xlabel(xlabels{d})
    if ~verticalAlignment
        ylabel({'Posterior';'probability density'})
    else
        if d==2
            ylabel({'Posterior probability density'})
        end
    end
    yl = ylim;
    if max(f)>0.9*yl(2)
        yl(2) = yl(2)*1.1;
        ylim(yl);
    end
    xl = xlim;
    xlim([min(0,xl(1)) xl(2)])
    if reverseLeave
        if d==1 || d==3
            xlim([xl(1) max(xl(2),0)])
        end
    end
    xl = xlim;
    if xl(1)
        line([0 0], yl, 'Color','k')
    end
    meanParticle = weights_all_linear' * thetas_all_linear(:,d); %signflip*meanThetas(d,winningModel(c),c,p);
    hdi95 = hdi95_weighted(thetas_all_linear(:,d),weights_all_linear);
    lower = meanParticle - hdi95(1); %signflip*hdi95s(d,1,winningModel(c),c,p);
    upper = hdi95(2) - meanParticle; %signflip*hdi95s(d,2,winningModel(c),c,p)
    errorbar(meanParticle, yl(2)*0.95, ...
        lower, upper, ...
        'horizontal', 'o','Color','k','MarkerFaceColor','k','MarkerSize',3)

    if forPub

        set(gca,'TickLength',[.03 .03],'TickDir','out')
        if verticalAlignment
            set(gca,'units','centimeters','position',[1,1+1.6*(4-d),2,0.6])
        else
            axis square
            set(gca,'units','centimeters','position',[3*(d-1)+1,1,2,2])
        end
        fontsize(6,'points')
        set(gca,'color','none')
        box off
    end
end


%% save naive params
save('naiveCombinedParams_theta_w_only.mat',...
    'thetas_all_linear','weights_all_linear');