
function [Rd, Rp] = calcLearningDepPot_abc(x,y,z,ta,tm,a)

mu_enter_0 = 1.4436 - 0.2885;


AvoidD = 1 + (1-x*(1-z))*y;
ApproachD = (1-x*(1-z))*(1+y);
AvoidD = AvoidD*a;
ApproachD = ApproachD*a;
Rd = normcdf( mu_enter_0 - (AvoidD-ApproachD)/sqrt( (ta + tm*AvoidD)^2 + (ta + tm*ApproachD)^2 ) ) ;

AvoidP = (z+(1-z)*x)*(1+y);
ApproachP = z + (z+(1-z)*x)*y;
AvoidP = AvoidP*a;
ApproachP = ApproachP*a;
Rp = normcdf( mu_enter_0 - (AvoidP-ApproachP)/sqrt( (ta + tm*AvoidP)^2 + (ta + tm*ApproachP)^2 ) ) ;

Rd = 1 - Rd;
Rp = 1 - Rp;
end