function t0 = u0_solve(ux0,vx0)
%U0_SOLVE
%    T0 = U0_SOLVE(UX0,VX0)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    27-Aug-2020 22:10:06

t2 = ux0.^2;
t3 = vx0.^2;
t4 = vx0.*1i;
t5 = -t2;
t6 = -t3;
t7 = -t4;
t9 = -1.0./(t4-ux0);
t8 = t7+ux0;
t10 = t5+t6;
t11 = sqrt(t10);
t12 = t9.*t11;
t0 = [log(t12).*-1i;log(t11./(t4-ux0)).*-1i];
