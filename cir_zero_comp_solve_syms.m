% File: cir_zero_comp_solve_syms.m

clearvars();

syms r c0 u0 v0 theta

comp0 = symfun( c0 + r*u0*cos(theta) + r*v0*sin(theta), [r c0 u0 v0 theta] );

eq0 = comp0( r, c0, u0, v0, theta );

soln = solve( 0==eq0, theta );

matlabFunction( soln, 'file', 'cir_zero_comp_solve.m' )

