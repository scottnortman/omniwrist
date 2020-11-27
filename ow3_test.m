% File:  ow3_test.m

clc();
close all;
clearvars();

% Create ow ow instance
%  ow = ow3( 'radius', 50, 'edge',40, 'alpha', pi/4, 'ang_len', 10);
 
  ow = ow3( 'radius', 40, 'edge',25, 'alpha', pi/2, 'ang_len',0);
%   ow = ow3( 'radius', 40, 'edge',25, 'alpha', pi/3, 'ang_len',0);
 
 tyler = 0.75;
 
 ow.plot3d('zoom', tyler);
 
rs = 2.75;
ts = 0;
 
 
 % Import 3d mouse class
 import mouse3D.*
 
 % Create a copy of the instance; this is needed otherwise weird stuffs
 hDrv = mouse3Ddrv;
 
 ow.drive( hDrv )
 
 
 %  while 1
%      
%      pause(0.01);
%      
%      [dec, az] = ow_spacemouse_input( hDrv, rs, ts );
%      
%      ow.drive( dec, az );
%      
%  end
%  

 
 
 
 
 
%  % Capture histograms 
%  tout = 5; % capture time, s
%  
%  tstart = tic();
%  tdelta = toc( tstart );
%  while tdelta < tout
%     
%      tdelta = toc( tstart )
%      
%      
%      
%  end

