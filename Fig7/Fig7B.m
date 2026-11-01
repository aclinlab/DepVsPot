load('naiveParticles.mat'); % dilutions 1-7
vinParams = naiveParticles(6).weights' * naiveParticles(6).theta; % dilution 6, i.e. 1e-1 vinegar
mu_enter = vinParams(1);
sigma_enter = vinParams(2);
mu_stay = vinParams(3);
sigma_stay = vinParams(4);


%% X_enter and X_stay
Xenter_samples = mu_enter+randn(1e6,1)*sigma_enter;
Xstay_samples = mu_stay+randn(1e6,1)*sigma_stay;

figure,hold on
[nxenter,xxenter] = hist(Xenter_samples,100);
[nxstay,xxstay] = hist(Xstay_samples,100);
plot(xxenter,nxenter/1e6)
plot(xxstay,nxstay/1e6)
plotcmap = lines(2);
line(mu_enter*[1 1], [0 max(nxenter/1e6)],'Color',plotcmap(1,:))
line(mu_enter*[1 1] + [0 sigma_enter], max(nxenter)/1e6*normpdf(1)/normpdf(0)*ones(2,1),'Color',plotcmap(1,:))
line(mu_stay*[1 1], [0 max(nxstay/1e6)],'Color',plotcmap(2,:))
line(mu_stay*[1 1] - [0 sigma_stay], max(nxstay)/1e6*normpdf(1)/normpdf(0)*ones(2,1),'Color',plotcmap(2,:))
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off
xlim([-4 4])



%% P(enter) and P(stay)
Penter_samples = normcdf(Xenter_samples);
[nenter,xenter] = hist(Penter_samples,100);
Pstay_samples = normcdf(Xstay_samples);
[nstay,xstay] = hist(Pstay_samples,100);
figure,hold on
plot(xenter,nenter/1e6)
plot(xstay,nstay/1e6)
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
axis square
fontsize(6,'points')
box off
