function [p1, p2] = cir_cir_int( cen1, r1, cen2, r2 )

p1 = [];
p2 = [];

nrm = norm( cen1 - cen2 );

if nrm > r1 + r2
   return; 
end

if nrm < abs(r1-r2)
   return; 
end

cosAlpha = (r1^2+nrm^2-r2^2)/(2*r1*nrm);

uAB = (cen2 - cen1)/nrm;
puAB = [uAB(2),-uAB(1)]';

p1 = cen1 + uAB*(r1*cosAlpha) + puAB*(r1*sqrt(1-cosAlpha^2));
p2 = cen1 + uAB*(r1*cosAlpha) - puAB*(r1*sqrt(1-cosAlpha^2));

end