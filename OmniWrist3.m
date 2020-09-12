classdef OmniWrist3 < handle
% CLASSDEF OmniWrist3 parallel mechanism class with four linkage arms;


	properties % TODO:  set access control
        
        % Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NUM_SIDES = 4;
        NUM_TXS = 6;
        
        TX0 = 1;    % Side base frame
        TX1 = 2;    % Base link frame
        TX2 = 3;    % Angle link frame
        TX3 = 4;    % Platform link frame
        TX4 = 5;    % Platform link base frame
        TX5 = 6;    % Platform frame
        
		% Main link parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mainLinkLength = 0;
		mainLinkRadius = 0;

		% Angle link parameters
		angleLinkLength = 0;
		angleLinkAngle = pi/4; %TODO remvoe hard coded value

		% Transforms
		linkTxs = [];
        
        % Joint angles
        qActual = [];
        qNominal = [];
        
        % Derived values
        apexLength = [];
        workspaceRadius = [];

        % Face, verticie data from STL for patch object 
		baseTx1FV 		= [];
		baseLinkTx2FV 	= [];
		angleLinkTx3FV 	= [];
		platLinkTx4FV 	= [];
		platformTx5FV 	= [];

		name = '';
        
        rb_tx_plat = eye(4,4);
		

	end % properties

	properties (GetAccess = protected, SetAccess = protected)

		world_hgt_rb = hgtransform();
		rb_hgt_tx = hgtransform();
        
        rb_tx_plat_hgt = hgtransform();


		
	end



	methods

		% Class constructor
		function obj = OmniWrist3( ...
			mainLinkLength,... 		% arg 1
			mainLinkRadius,...		% arg 2
			angleLinkLength,...		% arg 3
			angleLinkAngle, ...		% arg 4
			varargin )

			% Check input arguments (scalars )

			% Parse through variable argument list
			if ~isempty( varargin )

				for argnum=1:2:length(varargin)

					% Check for base STL
					if 		strcmpi( 'baseTx1STLFilePath', varargin{argnum} )
						obj.baseTx1FV = OmniWrist3.checkSTLFileAndLoad( varargin{argnum+1} );
					elseif 	strcmpi( 'baseLinkTx2STLFilePath', varargin{argnum} )
						obj.baseLinkTx2FV = OmniWrist3.checkSTLFileAndLoad( varargin{argnum+1} );
					elseif 	strcmpi( 'angleLinkTx3STLFilePath', varargin{argnum} )
						obj.angleLinkTx3FV = OmniWrist3.checkSTLFileAndLoad( varargin{argnum+1} );
					elseif 	strcmpi( 'platLinkTx4STLFilePath', varargin{argnum} )
						obj.platLinkTx4FV = OmniWrist3.checkSTLFileAndLoad( varargin{argnum+1} );
					elseif 	strcmpi( 'platformTx5STLFilePath', varargin{argnum} )
						obj.platformTx5FV = OmniWrist3.checkSTLFileAndLoad( varargin{argnum+1} );
					else
						warning( 'Unknown option passed...\n' );
					end
					
				end % for argnum=1:2:length(varargin)

			end % ~isempty( varargin )

			% If input args OK; assign to class members
			obj.mainLinkLength = mainLinkLength;
			obj.mainLinkRadius = mainLinkRadius;
			obj.angleLinkLength = angleLinkLength;
			obj.angleLinkAngle = angleLinkAngle;
            obj.apexLength = obj.mainLinkLength + obj.mainLinkRadius + (obj.angleLinkLength / tan(obj.angleLinkAngle/2));
            obj.workspaceRadius = 2*obj.apexLength*sin( obj.angleLinkAngle/2 );
            
            obj.qActual = zeros(4, obj.NUM_SIDES);
            obj.qNominal = zeros(4, obj.NUM_SIDES);
            
            obj.linkTxs = zeros(4,4,obj.NUM_TXS,obj.NUM_SIDES); %row,col,link,side;
            
            % TODO:  clean up initialization of frames.....
            qnom = [    pi/8,pi/2,-pi/2,-pi/8;
                        pi/8,pi/2,-pi/2,-pi/8;
                        pi/8,pi/2,-pi/2,-pi/8;
                        pi/8,pi/2,-pi/2,-pi/8; ]';

            obj.calculateTxs(qnom);

            % Calculate joint values given nominal pose; if a valid soln is
            % found internal parameters are updated per soln
           % obj.ikine( 0, 0 );
            

            obj.qNominal = qnom;

		end % OmniWrist3
        
        function qs = ikine( obj, dec, az )
        % METHOD ikine( obj, dec, az )
        %   
        % General approach:  The OmniWrist 3 mechanism is symmetric about a plane
        %	passing through all apex points for every given configuration.  This
        %	plane also always passes through the mid point of a vector from the mechanism
        %	origin (robot base frame 'rb') to the end effector platform origin. We
        %	can therefore solve for most of the kinematic joint angles given this
        %	fact by solving for plane-circle intersections.
        %
        %   tx = calc_ow3_tx_sdn( dec, az, r )
        %
        
        % Local variables
        qs = zeros( 4, obj.NUM_SIDES ); % save results to a local variable; then copy to obj

        % Given declination and azimuth, calculate corresponding 4x4
        % transform
        obj.rb_tx_plat = calc_ow3_tx_sdn( dec, az, obj.workspaceRadius );
       
    
       % Calc mid point and normal from 4x4 tx
       rb_uvec_mid = obj.rb_tx_plat(1:3,4) / norm( obj.rb_tx_plat(1:3,4) );
       rb_pt_mid = rb_uvec_mid * (obj.workspaceRadius/2);
    
        
        % given mid point and unit vector, project into each side TX0; note that TX0 frames are
        %	static w.r.t the base frame
        for side=1:obj.NUM_SIDES
            

 			%%% Solve for theta 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Solve wrt TX0; tx0_pt_mid = tx0_tx_rb * rb_pt_mid
            % Get global midpoint and plane normal uvec in local frame TX0
            tmp = obj.linkTxs(:,:,obj.TX0,side) \ [rb_pt_mid;1];
            tx0_pt_mid = tmp(1:3,1);
            tmp = obj.linkTxs(:,:,obj.TX0,side) \ [rb_uvec_mid;1];
            tx0_uvec= tmp(1:3,1);
            % Now we have plane defined w.r.t TX0 as tx0_pt_mid, tx0_uvec
            
            qs(1,side) = atan2( tx0_uvec(2), tx0_uvec(1) ) - pi/2 + (obj.angleLinkAngle/2);
            
            tx0_pt_apex = [obj.apexLength*cos( qs(1,side) ); obj.apexLength*sin( qs(1,side) ); 0];
            
% %             
% %             % Find plane - circle intersection point in local frame given
% %             % the following:
% %             %   Plane defined by point, normal: tx0_pt_mid, tx0_uvec
% %             %   Circle in XY plane defined by side N TX0 centered about the
% %             %   Z axis, with a radius of apexLength.  See 
% %             %   plane_circle_int_syms.m
% %             %   plane_circ_int.m
% %             %	[xp,yp,zp]' = point on plane [x,y,z]
% %             %	[xn,yn,zn]' = normal to plane [x,y,z]
% %             xp = tx0_pt_mid(1);
% %             yp = tx0_pt_mid(2);
% %             zp = tx0_pt_mid(3);
% %             xn = tx0_uvec(1);
% %             yn = tx0_uvec(2);
% %             zn = tx0_uvec(3);
% %             r = obj.apexLength;
% %             
% % % %             %%% BUG fix
% % % %             if abs(xn) < 1e-12
% % % %                 if xn < 0
% % % %                     xn = -1e-12;
% % % %                 else
% % % %                     xn=1e-12;
% % % %                 end
% % % %             end
% % % %             if abs(yn) < 1e-12
% % % %                 if yn < 0
% % % %                     yn = -1e-12;
% % % %                 else
% % % %                     yn=1e-12;
% % % %                 end
% % % %             end
% % % %             if abs(zn) < 1e-12
% % % %                 if zn < 0
% % % %                     zn = -1e-12;
% % % %                 else
% % % %                     zn=1e-12;
% % % %                 end
% % % %             end
% % % %             if abs(xp) < 1e-12
% % % %                 if xp < 0
% % % %                     xp = -1e-12;
% % % %                 else
% % % %                     xp=1e-12;
% % % %                 end
% % % %             end
% % % %             if abs(yp) < 1e-12
% % % %                 if yp < 0
% % % %                     yp = -1e-12;
% % % %                 else
% % % %                     yp=1e-12;
% % % %                 end
% % % %             end
% % % %             if abs(zp) < 1e-12
% % % %                 if zp < 0
% % % %                     zp = -1e-12;
% % % %                 else
% % % %                     zp=1e-12;
% % % %                 end
% % % %             end
% % % %             %%%%
% %             tx0_pts_soln = plane_circ_int( r, xn, xp, yn, yp, zn, zp );
% %             
% %             
% %        		
% %        		% Given circle plane intersection, there are two possible solutions; calculate angles for 
% %        		%	both; compare each against allowed ROM for theta 1 joint [-pi/2,pi/2]
% %             try
% %                 a1 = atan2( tx0_pts_soln(2,1), tx0_pts_soln(1,1) );
% %                 a2 = atan2( tx0_pts_soln(2,2), tx0_pts_soln(1,2) );            
% %                 % given two solns; pick one in range of [-pi/2,pi/2], store to
% %                 % local variable
% %                 if a1 > -pi/2 && a1 < pi/2
% %                     %obj.qActual( 1,side) = a1;
% %                     qs(1,side) = a1;
% %                     tx0_pt_apex = tx0_pts_soln(:,1);
% %                 else
% %                     qs(1,side) = a2;
% %                     %obj.qActual( 1,side) = a2;
% %                     tx0_pt_apex = tx0_pts_soln(:,2);
% %                 end
% %             catch
% %                 keyboard
% %             end
            
            

% %             
% %             %%% Solve for theta 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %             %Assumed that TX4 origin is coincident with the end effector platform
% %             %	origin.  For each side, assume that the corresponding link joint 
% %             %	frame is initially oriented the same as the TX0 frame (static) plus
% %             %	plus an additional rotation about the moving frame (TX0) Y axis by pi. 	
% %             % Solve for pt_apex wrt TX4; we need tx4_tx_t0
% %             %   tx4_pt_apex = tx4_tx_tx0 * tx0_pt_apex
% %             %tx4_pt_apex = inv(obj.linkTxs(:,:,obj.TX4,side )) * obj.linkTxs(:,:,obj.TX0,side) * [tx0_pt_apex;1];
% %             tx4_pt_apex = (obj.linkTxs(:,:,obj.TX4,side ) \ obj.linkTxs(:,:,obj.TX0,side)) * [tx0_pt_apex;1];
% %             qs(4,side) = atan2( tx4_pt_apex(2), tx4_pt_apex(1) )
% %             %qs(4,side) = - atan( tx4_pt_apex(2,1)/ tx4_pt_apex(1,1) );
% %             
            % Solve for theta 2
            % Project rb_pt_mid into tx1 frame
            %   tx1_pt_mid  = tx1_tx_rb * rb_pt_mid
            %               = inv( rb_tx_tx1 ) * rb_pt_mid
            
            
            % theta 4; project apex point into platform frame rotated by
            % rx((side-1)*(pi/2))
            plat_tx_tmp = obj.rb_tx_plat;
            plat_tx_tmp(1:3,1:3) = plat_tx_tmp(1:3,1:3)*OmniWrist3.rotz((side-1)*(pi/2));
            tmp_pt = (plat_tx_tmp \ obj.linkTxs(:,:,obj.TX0,side))*[tx0_pt_apex;1] ;
            qs(4,side) = atan2( tmp_pt(3), tmp_pt(1) );
            qs(4,side) = - qs(1,side);
            
            %tx1_pt_mid = inv( obj.linkTxs(:,:,obj.TX1,side) ) * [rb_pt_mid;1];
            tx1_pt_mid = obj.linkTxs(:,:,obj.TX1,side) \ [rb_pt_mid;1];
            % After projecting, theta2 is determined by the angle on the YZ
            % plane (ie X=0)
            qs(2,side) = atan2( tx1_pt_mid(3,1), tx1_pt_mid(2,1) ) + pi/2;
            
%           temp2 = acos( dot( obj.linkTxs(1:3,2,obj.TX1,side), obj.linkTxs(1:3,2,obj.TX2,side) ) )
            
            
            % Solve for theta 3
            qs(3,side) = qs(2,side) - pi;

        end 
        
        % Check possible soln against joint limits; if all are OK, assign
        % to object variable
        % TODO implement check
        obj.qActual = qs;
        
        % Re-calculate mechanism transforms given new valid angle set
        obj.calculateTxs( qs );
            
        end % ikine

		function calculateTxs( obj, qs )
        % METHOD calculateTxs( obj, qs )
        %   Calculates the mechanism 4x4 homogeneous transforms given a set
        %   of joint angles.

			% Calc frames; thse are w.r.t Tx0 frame
			for side = 1:obj.NUM_SIDES

				% Tx 0, rotate about X by pi/2, rb_tx_tx0
				obj.linkTxs(:,:,obj.TX0,side) = eye(4,4);
				obj.linkTxs(1:3,1:3,obj.TX0,side) = OmniWrist3.rotx(pi/2)*OmniWrist3.roty((side-1)*(pi/2));

				% q1, Tx 1, now rotate frame Tx1 about Z by q1, rb_tx_tx1
				obj.linkTxs(:,:,obj.TX1,side) = eye(4,4);
				obj.linkTxs(1:3,1:3,obj.TX1,side) = obj.linkTxs(1:3,1:3,obj.TX0,side)*OmniWrist3.rotz(qs(1,side));

				% q2, Tx 2, rotate about y by pi/2 and about z by q2, rb_tx_tx2
				obj.linkTxs(:,:,obj.TX2,side) = eye(4,4);
				obj.linkTxs(1:3,1:3,obj.TX2,side) = obj.linkTxs(1:3,1:3,obj.TX1,side)*OmniWrist3.roty(pi/2)*OmniWrist3.rotz(qs(2,side));

				% q3, Tx 3, rotate about y by angleLinkAngle then about z by q3,
				% translate per angle, rb_tx_tx3
				obj.linkTxs(:,:,obj.TX3,side) = eye(4,4);
				obj.linkTxs(1:3,1:3,obj.TX3,side) = obj.linkTxs(1:3,1:3,obj.TX2,side)*OmniWrist3.roty(-obj.angleLinkAngle)*OmniWrist3.rotz(qs(3,side));
				
				tx3_p_apex = [0 0 obj.apexLength]';
				tmp = obj.linkTxs(:,:,obj.TX2,side)*[tx3_p_apex;1];
				obj.linkTxs(1:3,4,obj.TX3,side) = tmp(1:3,1);            

				% q4, Tx 4, rb_tx_tx4, rotate about y by pi/2, then about z by q4
				obj.linkTxs(:,:,obj.TX4,side) = eye(4,4);
				obj.linkTxs(1:3,1:3,obj.TX4,side) = obj.linkTxs(1:3,1:3,obj.TX3,side)*OmniWrist3.roty(pi/2)*OmniWrist3.rotz(qs(4,side));
				tx4_p_apex = [0 0 -obj.apexLength]';
				tmp = obj.linkTxs(:,:,obj.TX3,side)*[tx4_p_apex;1];
				obj.linkTxs(1:3,4,obj.TX4,side) = tmp(1:3,1); 
                
                % TX5, rb_tx_tx5, rotate about x by -pi/2, then z by pi/2 *
                % side-1
                obj.linkTxs(:,:,obj.TX5,side) = obj.linkTxs(:,:,obj.TX4,side);
                %obj.linkTxs(1:3,1:3,obj.TX5,side) = obj.linkTxs(1:3,1:3,obj.TX4,side)*OmniWrist3.rotx(-pi/2)*OmniWrist3.rotz(pi+((side-1)*(pi/2)));
                obj.linkTxs(1:3,1:3,obj.TX5,side) = obj.linkTxs(1:3,1:3,obj.TX4,side)*OmniWrist3.rotz(pi)*OmniWrist3.rotx(pi/2)*OmniWrist3.rotz(-(side-1)*(pi/2));
                
			end % end for

		end % calculateTxs

		function axhan = plot3( obj, varargin )
		% METHOD axhan = plot3( obj, varargin )
		%
		%

			% Parameter defaults values
			defaultFrameScale = 10;	
            defaultFrameLineWidth = 2;	
            defaultKinematicLineWidth = 1;
            defaultPauseValue = 0;

            % Optional plot object default flags
            plotFramesFlag = 1;
            plotKineFlag = 0;
            plotTextFlag = 0;
            plotPathFlag = 0;
            plotStlFlag = 1;
            
            axhan = [];

            % Check number of passed arguments; note 'obj' is always passed
			if ~isempty( varargin )

				% Match 'name', 'value' pairs
				for argnum=1:2:length(varargin)
					if 		strcmpi('AxesHandle', varargin{argnum} )
						axhan = varargin{argnum+1};
					elseif 	strcmpi('FrameScale', varargin{argnum} )
						defaultFrameScale = varargin{argnum+1};
					elseif 	strcmpi('FrameLineWidth', varargin{argnum} )
						defaultFrameLineWidth = varargin{argnum+1};
					elseif 	strcmpi('KinematicLineWidth', varargin{argnum} )
						defaultKinematicLineWidth = varargin{argnum+1};
					elseif 	strcmpi('PauseValue', varagin{argnum} )
						defaultPauseValue = varargin{argnum+1};
					elseif strcmpi( 'PlotFramesFlag', varargin{argnum} )

					elseif strcmpi( 'PlotKineFlag', varargin{argnum} )

					elseif strcmpi( 'PlotTextFlag', varargin{argnum} )

					elseif strcmpi( 'PlotPathFlag', varargin{argnum} )
	
					elseif 	strcmpi( 'PlotSTLFlag', varargin{argnum} )
                        
                    elseif strcmpi( 'worldTxBase', varargin{argnum} )

					end	
				end

			end % ~isempty( varargin )

			if isempty( axhan )

				fhan = figure( 'name', obj.name, 'numbertitle', 'off' );
				axhan = axes( 'parent', fhan );
				axhan.NextPlot = 'replaceChildren';

				% Robot base tx wrt world
				obj.world_hgt_rb	= hgtransform( 'parent' , axhan, 'tag', 'world_hgt_rb');
                
                obj.rb_tx_plat_hgt = hgtransform( 'parent', obj.world_hgt_rb );
                
                obj.rb_tx_plat_hgt.Matrix = obj.rb_tx_plat;


				% Serial link 1 txs
				for side = 1:obj.NUM_SIDES
					for tx = 1:obj.NUM_TXS
						obj.rb_hgt_tx(tx, side) = hgtransform( 'parent', obj.world_hgt_rb );
					end
	            end

	        end % if isempty( axhan )
            
            % Update HGT matricies
            
            for side = 1:obj.NUM_SIDES
                for tx = 1:obj.NUM_TXS
                    obj.rb_hgt_tx(tx, side).Matrix = obj.linkTxs(:,:,tx,side);
                end
	        end
            

	        % Check flags to determine what should be plotted
	        if plotFramesFlag
                
     
                
                % Plot robot
                plot3( obj.rb_tx_plat_hgt, 0,0,0,'k*',...
                	[0 10], [0 0], [0 0],'r',...
                    [0 0], [0 10],[0 0], 'g',...
                    [0 0], [0 0], [0 10], 'b',...
                    'linewidth', 5,...
                    'tag', 'frame' );
% %                 
                %Plot robot base frame
%                 plot3(  obj.world_hgt_rb, 0,0,0,'k*',...
%                     [0 defaultFrameScale], [0 0], [0 0],'r',...
%                     [0 0], [0 defaultFrameScale],[0 0], 'g',...
%                     [0 0], [0 0], [0 defaultFrameScale], 'b',...
%                     'linewidth', defaultFrameLineWidth,...
%                     'tag', 'frame' );

	        	for side = 1:obj.NUM_SIDES
					for tx = 4%1:obj.NUM_TXS
	                    plot3(  obj.rb_hgt_tx(tx,side), 0,0,0,'k*',...
	                        [0 defaultFrameScale], [0 0], [0 0],'r',...
	                        [0 0], [0 defaultFrameScale],[0 0], 'g',...
	                        [0 0], [0 0], [0 defaultFrameScale], 'b',...
	                        'linewidth', defaultFrameLineWidth,...
	                        'tag', 'frame' );
                        hold on;
                    end
	            end
	        end % plotFramesFlag

	        if plotStlFlag
                
                if ~isempty( obj.baseTx1FV )
                    patch( 'parent', obj.rb_hgt_tx(1,side), obj.baseTx1FV, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
                end
                if ~isempty( obj.platformTx5FV )
                    patch( 'parent', obj.rb_hgt_tx(5,side), obj.platformTx5FV,'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
                end
                
	        	for side = [1,3]%1:obj.NUM_SIDES
                    if ~isempty( obj.baseLinkTx2FV )
                        patch( 'parent', obj.rb_hgt_tx(2,side), obj.baseLinkTx2FV, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
                    end
                    if ~isempty( obj.baseLinkTx2FV )
                        patch( 'parent', obj.rb_hgt_tx(3,side), obj.angleLinkTx3FV,'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
                    end
                    if ~isempty( obj.platLinkTx4FV )
                        patch( 'parent', obj.rb_hgt_tx(4,side), obj.platLinkTx4FV, 'facecolor', [1 1 1], 'edgecolor', [0 0 0], 'facealpha', 0.5 );
                    end
                end
            end
            
	       	axis equal vis3d;
            
            xlim([-150 150]);
            ylim([-150 150]);
            zlim([-150 150]);
            
            grid on

		end % plot3

	end %methods

	methods(Static)

		function r = rotx( ang )
			ca = cos( ang );
			sa = sin( ang );
			r = [   1   0   0   ;...
			        0   ca  -sa ;...
			        0   sa   ca ];
		end %rotx

		function r = roty( ang )
			ca = cos( ang );
			sa = sin( ang );	
			r = [	ca   0   sa;...
     				0    1   0;...
     				-sa  0   ca];
		end %roty

		function r = rotz( ang )
			ca = cos( ang );
			sa = sin( ang );
			r = [	ca   -sa   0;...
				    sa    ca   0;...
				    0     0    1];
		end %rotz

		function ret = checkSTLFileAndLoad( fpath )

			ret = [];

			if  2 ~= exist( fpath, 'file')
				warning( 'Invalid file [ %s ], not loading STL model...\n', fpath )
			else
				fprintf( 'Found file, attempting to load STL model [ %s ]...\n', fpath );
				try
					ret = stlread( fpath );
				catch
					warning('Unable to load file %s...\n', fpath );
				end
				if ~isempty(ret)
					fprintf('Loaded STL file %s...\n', fpath );
				end
			end

		end

	end %methods (Static)


end % classdef OmniWrist3