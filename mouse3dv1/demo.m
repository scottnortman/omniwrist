 close all; clc
 clear

import mouse3D.*


hDrv = mouse3Ddrv; %instantiate driver object
hMon = mouse3Dmonitor(hDrv); %use of the monitor is optional
%% Make a 3D image
 hFig = figure('Color', [0 0 0]);
 hAxes = axes;
 surf(hAxes,peaks); 
 shading interp
 colormap bone
 axis off; 

%% Make a camera control object for the figure
 j = mouse3DfigCtrl(hDrv,hAxes);
