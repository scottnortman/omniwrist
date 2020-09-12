% File:  OmniWrist3Syms.m

clear all;

syms 	mainLinkLength mainLinkRadius angleLinkLength angleLinkAngle sideAngle...
		... Frame joint angles
		side1theta1 side1theta2 side1theta3 side1theta4...
		side2theta1 side2theta2 side2theta3 side2theta4...
		side3theta1 side3theta2 side3theta3 side3theta4...
		side4theta1 side4theta2 side4theta3 side4theta4...
		...
		apexLength rx ry rz r t temp...
        ...
        side1tx0 side1tx1 side1tx2 side1tx3 side1tx4...



rotx2tr = symfun( [1 0 0 0; 0 cos(rx) -sin(rx) 0; 0 sin(rx) cos(rx) 0; 0 0 0 1], rx );
roty2tr = symfun( [cos(ry) 0 sin(ry) 0; 0 1 0 0; -sin(ry) 0 cos(ry) 0; 0 0 0 1], ry );
rotz2tr = symfun( [cos(rz) -sin(rz) 0 0; sin(rz) cos(rz) 0 0; 0 0 1 0; 0 0 0 1], rz );

% apexLength = mainLinkLength + mainLinkRadius + (angleLinkLength/tan(angleLinkAngle/2));
% 
% 
% side = 1;
% side1tx0 = rotx2tr(pi/2)*roty2tr((side-1)*pi/2);
% side1tx1 = side1tx0*rotz2tr(side1theta1);
% side1tx2 = side1tx1*roty2tr(pi/2)*rotz2tr(side1theta2);
% side1tx3 = side1tx2*roty2tr(-angleLinkAngle)*rotz2tr(side1theta3);
% temp = side1tx3*[0 0 -apexLength 1]'; 


ori = 


