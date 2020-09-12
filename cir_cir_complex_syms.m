% file:  circle_circle_intersection_3d_syms.m

%Given the definition of a parametric circle in 3d space with each circle
%defined by:
%
%       p = c + u.*r.*cos(theta) + v.*r.*sin(theta)
%
%   where
%
%       p, 3x1 vector point on circle
%       c, 3x1 center of circle
%       u, 3x1 unit vector in the plane of circle, perpendicular to axis
%       v, 3x1 unit vector in the plane of circle, perpendicular to u
%       theta, [0,2*pi], variable parameter

clearvars();

syms    px0 py0 pz0 ...
        cx0 cy0 cz0 ...
        ux0 uy0 uz0 ...
        vx0 vy0 vz0 ...
        theta0 ...
        r0 ...
        p0...
        c0...
        u0...
        v0...
        ...
        px1 py1 pz1 ...
        cx1 cy1 cz1 ...
        ux1 uy1 uz1 ...
        vx1 vy1 vz1 ...
        theta1 ...
        r1 ...
        p1...
        c1...
        u1...
        v1
   
    p0 = [px0;py0;pz0];
    c0 = [cx0;cy0;cz0];
    u0 = [ux0;uy0;uz0];
    v0 = [vx0;vy0;vz0];
    
    p1 = [px1;py1;pz1];
    c1 = [cx1;cy1;cz1];
    u1 = [ux1;uy1;uz1];
    v1 = [vx1;vy1;vz1];
    
    
    
    circle0 = symfun( p0 == c0 + u0*r0*cos(theta0) + v0*r0*sin(theta0), ...
        [px0 py0 pz0 cx0 cy0 cz0 ux0 uy0 uz0 vx0 vy0 vz0 r0 theta0] );
    
    circle1 = symfun( p1 == c1 +  u1*r1*cos(theta1) + v1*r1*sin(theta1), ...
        [px1 py1 pz1 cx1 cy1 cz1 ux1 uy1 uz1 vx1 vy1 vz1 r1 theta1] );
    
    eq0 = circle0( px0, py0, pz0, cx0, cy0, cz0, ux0, uy0, uz0, vx0, vy0, vz0, r0, theta0 )
    eq1 = circle1( px1, py1, pz1, cx1, cy1, cz1, ux1, uy1, uz1, vx1, vy1, vz1, r1, theta1 )
    
    soln = solve( eq0==eq1, [theta0, theta1] )
    
    t0 = solve( 0==r0*ux0*cos(theta0) + r0*vx0*sin(theta0), theta0 );
    
    matlabFunction( t0, 'file', 'u0_solve.m' )
    
    
    
        






