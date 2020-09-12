% File:  ow3_test.m

close all;
clearvars();


% Mechanism parameters
linkLength = 45;
linkRadius = 40;
angleLinkLength = 20;
angleLinkAngle = pi/4;

% Define path to STL files
STLBaseDir = 'D:\Managed\MagicLeap\emtracking\omniwrist-z\';
baseTx1STLFilePath 		= [STLBaseDir, 'BLOCK_END_EFFECTOR_TX1.STL'];
baseLinkTx2STLFilePath 	= [STLBaseDir, 'BRACKET_CIRCLE_ASSY_TX2.STL'];
angleLinkTx3STLFilePath = [STLBaseDir, 'BRACKET_ANGLE_ASSY_TX3.STL'];
platLinkTx4STLFilePath 	= [STLBaseDir, 'BRACKET_CIRCLE_ASSY_TX4.STL'];
platformTx5STLFilePath 	= [STLBaseDir, 'BLOCK_END_EFFECTOR_TX5.STL'];

ow3 = OmniWrist3( linkLength, ...
	linkRadius, ...
	angleLinkLength, ...
	angleLinkAngle, ...
	'baseTx1STLFilePath', baseTx1STLFilePath, ...
	'baseLinkTx2STLFilePath', baseLinkTx2STLFilePath,...
	'angleLinkTx3STLFilePath', angleLinkTx3STLFilePath,...
	'platLinkTx4STLFilePath', platLinkTx4STLFilePath, ...
	'platformTx5STLFilePath', platformTx5STLFilePath...
    );


qnom = [    pi/8,pi/2,-pi/2,-pi/8;
            pi/8,pi/2,-pi/2,-pi/8;
            pi/8,pi/2,-pi/2,-pi/8;
            pi/8,pi/2,-pi/2,-pi/8; ]';
