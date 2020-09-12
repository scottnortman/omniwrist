% 	File:  	ow_kine_gui.m
%	Desc:	Basic GUI to model the base and platform kinematics of the OmnWrist mechanism
function [] = ow_kine_gui()

%% Initialization
clearvars();
close all;
clc();


%% File constants
AXS = 25;
LIM = 250;
RADIUS = 100;       % distance from base origin to platform origin
ALPHA= pi/4;   % assumed fixed for now




STLBaseDir = 'D:\Managed\MagicLeap\emtracking\omniwrist-z\';

%% Set up figure and axis
fhan = figure(  'name', 'OmniWristIII Kinematics', ...
                'numbertitle', 'off',...
                'MenuBar', 'none' );
            
axhan = axes( 'parent', fhan, ...
    'projection', 'perspective',...
	'position', [0.05 0.2 0.7 0.6] );
axhan.Toolbar.Visible = 'on';
set( axhan, 'nextplot', 'replacechildren' );
axis vis3d equal 
grid on
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim( [-LIM,LIM] );
ylim( [-LIM,LIM] );
zlim( [-LIM,LIM] );
hold on


% Create hgt for later updating of graphics
hg0 = hgtransform( 'parent', axhan );           % Base frame
hgPlatform = hgtransform( 'parent', axhan );    % Platform frame
hgMidpoint = hgtransform( 'parent', axhan );    % Midpoint frame for kinematics

% Calc. TX5 (end effector) frame given starting dec, az, r
hgPlatform.Matrix = calc_ow3_tx_sdn( 0, 0, RADIUS );

% Plot base frame (same a TX0)
plot3( hg0, ...
    [0 AXS], [0 0], [0 0], 'r',...
    [0 0], [0 AXS], [0 0], 'g',...
    [0 0], [0 0], [0 AXS], 'b', ...
    0,0,0,'k*', 'linewidth', 1 );

% plot initial end effector platform frame
plot3( hgPlatform, ...
    [0 AXS], [0 0], [0 0], 'r',...
    [0 0], [0 AXS], [0 0], 'g',...
    [0 0], [0 0], [0 AXS], 'b', ...
    0,0,0,'k*', 'linewidth', 3 );

% Given RADIUS, we have a right triangle with the 90° corner at the
% midpoint, a side from said point to the angle apex formed by the angled
% links, a side from said midpoint to the base origin, and the hypotenuse
% from the origin to said apex point.
%  Also note each aforementinoed leg triangle is upper / lower symmetric
%   about an over defined plane containing all apex points.
% Note that there are 4 apex points which can be referenced by the platform
% coordinate frame in the -X, +X, -Y, +Y directions.  The length of the
% line from the midpoint is
lenMidToApex = (RADIUS/2) / (tan(ALPHA/2));


% Plot midpoint frame
plot3( hgMidpoint, ...
    [0 AXS], [0 0], [0 0], 'r',...
    [0 0], [0 AXS], [0 0], 'g',...
    [0 0], [0 0], [0 AXS], 'b', ...
    0,0,0,'k*', 'linewidth', 2);

% Plot apex points and lines
plot3( hgMidpoint, ...
    [0 lenMidToApex], [0 0], [0 0], 'k--', ...
    [0 0], [0 lenMidToApex], [0 0], 'k--',...
    [0 -lenMidToApex], [0 0], [0 0], 'k--', ...
    [0 0], [0 -lenMidToApex], [0 0], 'k--',...
    lenMidToApex, 0, 0, 'k*',...
    0, lenMidToApex, 0, 'k*', ...
    -lenMidToApex, 0, 0, 'k*',...
    0, -lenMidToApex, 0, 'k*' );

    

%% Load STL for visualization
platformTxPlatFV = stlread( [STLBaseDir, 'BLOCK_END_EFFECTOR_TX_PLAT.STL'] );
patch( 'Parent', hgPlatform, 'Vertices', platformTxPlatFV.Points, 'Faces', platformTxPlatFV.ConnectivityList, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );


%% UI elements
pnl = uipanel( 'title', 'Control', ...
	'units', 'normalized',...
	'position', [0.775 0 0.2 1] );

reset_btn = uicontrol( pnl, ...
    'style', 'pushbutton', ...
    'units', 'normalized', ...
    'position', [0 0.95, 1, 0.05], ...
    'string', 'RESET', ...
    'callback', @resetCallback );

%Azimuth slider
slider_azimuth = uicontrol( pnl,...
    'style', 'slider', ...
    'units', 'normalized',...
    'min', -pi, 'max', pi,...
    'value', 0, ...
    'position', [0 0.025 1, 0.4],...
    'callback', @azimuthCallback );
    
%Declination slider
slider_declination = uicontrol( pnl,...
    'style', 'slider', ...
    'units', 'normalized',...
    'min', 0, 'max', pi/2,...
    'value', 0, ...
    'position', [0 0.5 1, 0.4],...
    'callback', @declinationCallback );

% Azimuth text value
azimuth_text = uicontrol( pnl, ...
    'style', 'text', ...
    'enable', 'inactive', ...
    'string', 'Azimuth' , ...
    'units', 'normalized', ...
    'fontweight', 'bold',...
    'position', [0 0.425, 1, 0.05]);
    

% Declination text value
declination_text = uicontrol( pnl, ...
    'style', 'text', ...
    'enable', 'inactive', ...
    'string', 'Declination' , ...
    'units', 'normalized', ...
    'fontweight', 'bold',...
    'position', [0 0.9, 1, 0.05] );

radLineHan = plot3(axhan,  [0 hgPlatform.Matrix(4,1)], [0 hgPlatform.Matrix(2,4)], [0 hgPlatform.Matrix(3,4)], 'm--', 'LineWidth', 2);
dd = hgPlatform.Matrix(1:3,4)./RADIUS; % unit vector from origin to center

azimuthCallback();
declinationCallback();





%% Callbacks

    function resetCallback( varargin )
        

        slider_azimuth.Value = 0;
        slider_declination.Value = 0;
        
        azimuthCallback();
        declinationCallback();
        
    end

    function azimuthCallback( varargin  )
        azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
        axes_update();
    end


    function declinationCallback( varargin )
        declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );
        axes_update();
        
    end

    function axes_update()
        [tx, txMid] = calc_ow3_tx_sdn( slider_declination.Value, ...
            slider_azimuth.Value, RADIUS );
        
        tx
        
        hgPlatform.Matrix = tx;
        
        
        %Plot line from base frame to end effector platform
        radLineHan.XData = [0 tx(1, 4)];
        radLineHan.YData = [0 tx(2, 4)];
        radLineHan.ZData = [0 tx(3, 4)];


        hgMidpoint.Matrix = txMid;  

        
        
    end

    
end




