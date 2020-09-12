% file: plane_circle_int_syms.c
%

% Given a plane defined by a point [xp,yp,zp]' and a normal [xn,yn,zn]' we
% can write the equation as
% xn*(x-xp) + yn*(y-yp) + zn*(z-zp) = 0
% 
% We can assume that the intersection circle is in the XY plane (by design)
% so we let z = 0, therefore the above equation can simplify to:
%   xn*(x-xp) + yn*(y-yp) + zn*(-zp) = 0
%
%  We also have a circle centered at [xc,yc]' with a radius r, in the XY
%  plane therefore the equation is 
%   (x-xc)^2 + (y-yc)^2 = r^2
%
% Given the two above equations:
%   (1) xn*(x-xp) + yn*(y-yp) + zn*(-zp) = 0
%   (2) (x-xc)^2 + (y-yc)^2 = r^2
% we need to solve for the two unkowns x,y

clearvars();

syms xp yp zp xn yn zn xc yc r x y

plane_z0 = symfun( xn*(x-xp) + yn*(y-yp) + (zn*-zp) == 0, [xn yn zn xp yp zp x y] );
circle = symfun( (x-xc)^2 + (y-yc)^2 == r^2, [xc yc x y r] );
% To simplify, we know that xc = 0 and yc = 0 (ie centered at origin)
circle_c0 = symfun( x^2 + y^2 == r^2, [x y r] );

eq1 = plane_z0( xn, yn, zn, xp, yp, zp, x, y );
%eq2 = circle( xc, yc, x, y, r );
eq2 = circle_c0( x, y, r );

soln = solve( [eq1, eq2], [x, y] );

% Two solutions, each as column vectors
pta = [soln.x(1);soln.y(1)];
ptb = [soln.x(2);soln.y(2)];
pts = [pta,ptb];

% Convert solutions to callable MATLAB functions
matlabFunction( pts, 'file', 'plane_circ_int.m' );


