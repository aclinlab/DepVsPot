y = 7/3; %overlapFraction = 0.7 % y is the fraction of overlapping cells divided by the fraction of unique cells to A
z = 0; % depression goes from 1 to z; potentiation goes from z to 1
sValues = 0:0.01:1; % s is the "plasticity fraction"
[Rd, Rp] = deal(zeros(length(sValues),1));
ta = 0.0215;
tm = 0.171;
a = 0.06; % 0.2 cells respond to A, 0.06 unique to A, 0.14 also respond to B

for si=1:length(sValues)
    s = sValues(si);
    AvoidD = 1 + (1-s*(1-z))*y;
    ApproachD = (1-s*(1-z))*(1+y);
    AvoidD = AvoidD*a;
    ApproachD = ApproachD*a;
    % Ta + Tm*mu multiplicative noise
    Rd(si) = normcdf( (AvoidD-ApproachD)/sqrt( (ta + tm*AvoidD)^2 + (ta + tm*ApproachD)^2 ) ) ;

    AvoidP = (z+(1-z)*s)*(1+y);
    ApproachP = z + (z+(1-z)*s)*y;
    AvoidP = AvoidP*a;
    ApproachP = ApproachP*a;
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

