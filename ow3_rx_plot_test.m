% Script to test the ow3 rotation calculation

clearvars();
close all;
clc();

AXS = 25;
LIM = 250;
RADIUS = 100;

% Set up figure, homogeneous TX
fhan=figure();
ahan=axes('parent', fhan);
set( ahan, 'nextplot', 'replacechildren' );
axis vis3d equal
grid on
hg0 = hgtransform( 'parent', ahan );
hgstl = hgtransform( 'parent', ahan );

% Calc. TX5 (end effector) frame given starting dec, az, r
hgstl.Matrix = calc_ow3_tx_sdn( 0, 0, RADIUS );

% Plot base frame (same a TX0)
plot3( hg0, ...
    [0 AXS], [0 0], [0 0], 'r',...
    [0 0], [0 AXS], [0 0], 'g',...
    [0 0], [0 0], [0 AXS], 'b', ...
    0,0,0,'k*', 'linewidth', 1 );


plot3( hgstl, ...
    [0 AXS], [0 0], [0 0], 'r',...
    [0 0], [0 AXS], [0 0], 'g',...
    [0 0], [0 0], [0 AXS], 'b', ...
    0,0,0,'k*', 'linewidth', 3 );



xlim(LIM*[-1 1]);
ylim(LIM*[-1 1]);
zlim(LIM*[-1 1]);

STLBaseDir = 'D:\Managed\MagicLeap\emtracking\omniwrist-z\';
%platformTx5FV 	= stlread( [STLBaseDir, 'BLOCK_END_EFFECTOR_TX5.STL'] );
%patch( 'Parent', hgstl, platformTx5FV, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5 );
%patch( 'Parent', hgstl, 'Vertices', platformTx5FV.Points, 'Faces', platformTx5FV.ConnectivityList, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5 );
%frameFV = stlread( [STLBaseDir, 'frame.STL'] );
%patch( 'Parent', hgstl, 'Vertices', frameFV.Points, 'Faces', frameFV.ConnectivityList, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
platformTxPlatFV = stlread( [STLBaseDir, 'BLOCK_END_EFFECTOR_TX_PLAT.STL'] );
patch( 'Parent', hgstl, 'Vertices', platformTxPlatFV.Points, 'Faces', platformTxPlatFV.ConnectivityList, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );

NSTEPS = 1000;

%dec = pi/3*ones(NSTEPS);
dec =0.5*pi*sin(linspace(0, 10*pi, NSTEPS) );
%dec =0.5*pi*sin(linspace(0, 0, NSTEPS) );
az  = linspace( 1, 1, NSTEPS );
%az = linspace( 0, 0, NSTEPS );


% movie
mov(1:NSTEPS) = struct('cdata', [],...
                        'colormap', []);

input('Y?');


for step = 1:NSTEPS
   
    tx = calc_ow3_tx_sdn( dec(step), az( step),100 );
    %hgstl.Matrix = tx;
    %tmp = tx
    %tmp(1:3,1:3) = tmp(1:3,1:3)*rotx(pi/2);
    %hgstl.Matrix = tmp;
    hgstl.Matrix = tx;

%  xlim(LIM*[-1 1]);
% ylim(LIM*[-1 1]);
% zlim(LIM*[-1 1])

   %mov(step) = getframe(ahan);
   %pause(0.025);
    
   drawnow();
end

% movie2avi(mov, 'ow3.avi', 'compression', 'None');


