% For equal starting points: depression goes from 1 to z; potentiation goes from 1 to 2-z. 
% For equal dynamic ranges: depression goes from 1 to z; potentiation goes from z to 1. 
z = 0; 
sValues = 0:0.01:1; % s is the "plasticity fraction"
ta = 0.0215;
tm = 0.171;
a = 0.06; % 0.2 cells respond to A, 0.06 unique to A, 0.14 also respond to B. Also 0.06 are unique to B.
c = 0.14;
b = 0.06;

[AvoidD, ApproachD, AvoidP, ApproachP] = deal(zeros(length(sValues),3,2)); % s, single/discr/generalization, sameStart vs sameDynamicRange

startForPot = [1, z];
for i=1:3
    for j=1:2
        stPot = startForPot(j);
        for si=1:length(sValues)
            s=sValues(si);
            switch i
                case 1 % single odor
                    AvoidD(si,i,j) = a + c;
                    ApproachD(si,i,j) = (1-s*(1-z))*(a+c);
                    AvoidP(si,i,j) = (stPot+s*(1-z))*(a+c);
                    ApproachP(si,i,j) = stPot*(a+c);
                case 2 % A punished, B rewarded: test avoidance of A
                    AvoidD(si,i,j) = a + (1-s*(1-z))*c;
                    ApproachD(si,i,j) = (1-s*(1-z))*(a+c);
                    AvoidP(si,i,j) = (stPot+(1-z)*s)*(a+c);
                    ApproachP(si,i,j) = stPot*a + (stPot+(1-z)*s)*c;
                case 3 % A punished, B neutral: test generalization to B
                    AvoidD(si,i,j) = b + c;
                    ApproachD(si,i,j) = b + (1-s*(1-z))*c;
                    AvoidP(si,i,j) = stPot*b + (stPot+s*(1-z))*c;
                    ApproachP(si,i,j) = stPot*(b+c);
            end
        end
    end
end

Rd = normcdf( (AvoidD-ApproachD)./sqrt( (ta + tm.*AvoidD).^2 + (ta + tm.*ApproachD).^2 ) ) ;
Rp = normcdf( (AvoidP-ApproachP)./sqrt( (ta + tm.*AvoidP).^2 + (ta + tm.*ApproachP).^2 ) ) ;

ylabels = {'P(avoid)', 'P(avoid), Odor A', 'P(avoid), Odor B'};
figure,
for i=1:3
    for j=1:2
        subplot(3,2,j+(i-1)*2)
        hold on
        plot(sValues,Rd(:,i,j),'Color','r');
        plot(sValues,Rp(:,i,j),'Color','b');
        set(gca,'units','centimeters','position',[1+(j-1)*3,1+(4-i)*3,2,2])
        set(gca,'TickLength',[.03 .03],'TickDir','out')
        ylim([0.5 1]);
        yticks([0.5 0.75 1]);
        xlabel('Plasticity progress')
        ylabel(ylabels{i})
        axis square
        fontsize(6,'points')
        box off

    end
end
