z = 0; % depression goes from 1 to z; potentiation goes from z to 1
sValues = 0:0.01:1; % s is the "plasticity fraction"
[Rd, Rp] = deal(zeros(length(sValues),1));
ta = 0.0215; % T_a and T_m values from Fig 7
tm = 0.1710;
n = 0.2; % fraction of KCs activated by the odor

for si=1:length(sValues)
    s = sValues(si);
    AvoidD = n; 
    ApproachD = n*(1-s*(1-z)); 
    % Ta + Tm*mu multiplicative noise
    Rd(si) = normcdf( (AvoidD-ApproachD)/sqrt( (ta + tm*AvoidD)^2 + (ta + tm*ApproachD)^2 ) ) ;

    AvoidP = n*(z+(1-z)*s); %*(1+y);
    ApproachP = z*n;% + (z+(1-z)*s)*y;
    % Ta + Tm*mu multiplicative noise
    Rp(si) = normcdf( (AvoidP-ApproachP)/sqrt( (ta + tm*AvoidP)^2 + (ta + tm*ApproachP)^2 ) ) ;

end

figure, hold on
plot(sValues,Rd,'Color','r');
plot(sValues,Rp,'Color','b');
set(gca,'units','centimeters','position',[1,1,2,2])
set(gca,'TickLength',[.03 .03],'TickDir','out')
ylim([0.5 1]);
yticks([0.5 0.75 1]);
xlabel('Plasticity progress')
ylabel('P(avoid)')
axis square
fontsize(6,'points')
box off

