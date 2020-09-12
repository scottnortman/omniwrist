classdef ow3 < handle
    
    properties
        
        % Instance constants
        radius = [];        % Distance between kinematic centers
        edge = [];          % Edge length of base and platform squares
        alpha = [];         % radians, total angle between upper/lower arms
        ang_len =[];        % Length of edge of angle link (optional)
        name = 'ow3';
        
        % Homogeneous transforms of the actual instance state
        %world_tx_base = eye(4,4);   
        %base_tx_plat = eye(4,4);
        %base_tx_mid = eye(4,4);
        % Apex points of
        %base_pt_apex= eye(4,4);
        
        fighan = [];
        axhan = [];
       
        
        AXLEN = [];
        
        
        
        
    end
    
    methods
        
        function obj = ow3( varargin )
            % FUNCTION obj = ow3( varargin )
            %   Class constructor method for the OmniWrist3 object.
            %
            
            try
            
                for argnum=1:2:nargin

                    name = varargin{argnum};
                    value = varargin{argnum+1};

                    switch name

                        case 'radius'
                            obj.radius = value;
                        case 'edge'
                            obj.edge = value;
                        case 'alpha'
                            obj.alpha = value;
                        case 'name'
                            obj.name = value;
                        case 'ang_len'
                            obj.ang_len = value; 
                        otherwise
                            error('Invalid "name, value" pair');
       
                    end % switch

                end % for
                
            catch ME
                error( ME.identifier, 'ow3:ow3:args' );
                rethrow( ME );
            end
            
            
            %TODO: assert properties
            
            % 
            
            obj.AXLEN = 0.1*obj.radius;
            
            % if ang_len not passed, assume size
            if isempty( obj.ang_len )
                obj.ang_len = 0.25*obj.radius;
            end
            
            
        end %ow3
       
        function [q, ikine_soln] = ikine( obj, declination, azimuth )
            % FUNCTION [q, ikine_soln] = ikine( declination, azimuth )
            %   Inverse kinematic calculation for the ow3 mechanism.  Given
            %   the desired declination and azimuth pointing angles,
            %   calculates the corresponding joint angle solution.
            %
            %   INPUTS
            %       declination,    radians, [0, pi/2], 1 x n
            %       azimuth,        radians, [-pi, pi], 1 x n
            %       
            %   OUTPUTS
            %       q,              rad, 4 x n col vector, joints
            %       ikine_soln      array of structures
            %
            %   EXAMPLE
            %       >>  q = ow3.ikine( 0, 0 ); % Nominal zero
            %
            %   NOTE
            %       
            %       When in the nominal zero position, the orientation of
            %       the base frame and platform frame are identical, with
            %       the platform frame origin translated along the base
            %       frame Z axis a positive distance of 'radius'.
            %
            %       Joint 1 rotates about the positive base X axis, and
            %       joint 2 rotates about the positive base Y axis.  
            %
            %       The returned data strture list ikine_soln contains and
            %       array of datat structures, each with tansforms corresponding to
            %       the inverse kinematic solution for each input value
            %       set.  Included transforms are all with respect to the
            %       mechanism base frame.
            %
            %           base_tx_plat,       end effector platform
            %           base_tx_link1_bot,  joint 1 bottom link
            %           base_tx_link1_ang,  joint 1 angle link
            %           base_tx_link1_top,  joint 1 top link
            %           base_tx_link2_bot,  joint 1 bottom link
            %           base_tx_link2_ang,  joint 1 angle link
            %           base_tx_link2_top,  joint 1 top link
            %           base_tx_link3_bot,  joint 1 bottom link
            %           base_tx_link3_ang,  joint 1 angle link
            %           base_tx_link3_top,  joint 1 top link
            %           base_tx_link4_bot,  joint 1 bottom link
            %           base_tx_link4_ang,  joint 1 angle link
            %           base_tx_link4_top,  joint 1 top link
            %
            %       
            %       

                if length(declination) ~= length(azimuth)
                    error( 'Length of arguments much match.');
                end 
                
                
                % Define 4 apex points corresponding to the axis intersections 
                % of the angled linkages. These points are always in the
                % XY plane of the txMid coordinate frame, with  a fixed offset
                % from the txMid origin, determined by the angle alpha
                offset = (obj.radius/2) / tan(obj.alpha/2);
                
                mid_pt_apex =   [   offset,     0,          0;
                                    0,          offset,     0;
                                    -offset,    0,          0;
                                    0,          -offset,    0]';
                
                
                                
                len = length(declination);                
                q = zeros(4, len );
                ikine_soln = struct( ...
                    'q', cell(1,len),...
                    'base_pt_apex', cell(1,len),...
                    'base_tx_mid', cell(1,len),...
                    'base_tx_plat', cell(1,len),...
                    'base_tx_link1_bot', cell( 1, len ),...
                    'base_tx_link1_mid', cell( 1, len ),...
                    'base_tx_link1_top', cell( 1, len ),...
                    'base_tx_link2_bot', cell( 1, len ),...
                    'base_tx_link2_mid', cell( 1, len ),...
                    'base_tx_link2_top', cell( 1, len ),...
                    'base_tx_link3_bot', cell( 1, len ),...
                    'base_tx_link3_mid', cell( 1, len ),...
                    'base_tx_link3_top', cell( 1, len ),...
                    'base_tx_link4_bot', cell( 1, len ),...
                    'base_tx_link4_mid', cell( 1, len ),...
                    'base_tx_link4_top', cell( 1, len ) );
                                
                for ii=1:length( declination )
                    
                    dec = declination( ii );
                    az = azimuth( ii );
                    
                    % Given the desired declination angle and azimuth, the mechanism end
                    % effector platform origin translates along the sphere surface radius, 
                    % subtending half of the angle that the platform rotates.  
                    % Additionally Matlab requires the input as elevation.
                    [x,y,z] = sph2cart( az, (pi/2)-(dec/2), obj.radius );
                    
                    % Caclulate the Z axis orientation as a function of azimuth and
                    % declination; Note that this angle varies from 0 to pi/2 as the frame
                    % origin point subtends an angle of 0 to pi/4, ie half of the
                    % orientation angle, along the sphere surface per above
                    base_txLocal_plat=eye(4,4);
                    base_txLocal_plat(1:3,1:3) = rotz( az ) * roty( dec ) * rotz( -az );
                    base_txLocal_plat(1:3,4) = [x y z]';
                    
                    % calc mid frame
                    [xmid,ymid,zmid] = sph2cart( az, (pi/2)-(dec/2), obj.radius/2 );
                    base_txLocal_mid=eye(4,4);
                    base_txLocal_mid(1:3,1:3) = rotz( az ) * roty( dec/2 ) * rotz( -az );
                    base_txLocal_mid(1:3,4) = [xmid,ymid,zmid]';
                    ikine_soln(ii).base_tx_mid = base_txLocal_mid;
                    
                    % Transform above points into the base frame to solve for
                    % primary and secondary joint values wrt the base frame
                    base_ptLocal_apex = zeros(3,4);
                    for jj=1:4
                        pp = base_txLocal_mid * [mid_pt_apex(:,jj);1];
                        base_ptLocal_apex(:,jj) = pp(1:3,1);
                    end
                    %debug

                    % Apex points in the base frame
                    ikine_soln(ii).base_pt_apex = base_ptLocal_apex;
                    
                    % Create static base coordinate frames for each joint
                    base_tx_j1_fixed = eye(4,4);
                    base_tx_j1_fixed(1:3,1:3) = roty(pi/2) * rotz(pi/2);
                    base_tx_j1_fixed(1:3,4) = [obj.edge/2, 0, 0]';
                    
                    base_tx_j2_fixed = eye(4,4);
                    base_tx_j2_fixed(1:3,1:3) = rotx(-pi/2) * rotz(pi);
                    base_tx_j2_fixed(1:3,4) = [0 obj.edge/2, 0]';
                    
                    base_tx_j3_fixed = eye(4,4);
                    base_tx_j3_fixed(1:3,1:3) = roty(-pi/2) * rotz(-pi/2);
                    base_tx_j3_fixed(1:3,4) = [-obj.edge/2, 0, 0]';
                    
                    base_tx_j4_fixed = eye(4,4);
                    base_tx_j4_fixed(1:3,1:3) = rotx(pi/2);
                    base_tx_j4_fixed(1:3,4) = [0, -obj.edge/2, 0]';
                    
                    % Transform apex points into corresponding local joint
                    % frames
                    j1_static_pt_apex2 = base_tx_j1_fixed \ [base_ptLocal_apex(:,2);1];
                    j2_static_pt_apex3 = base_tx_j2_fixed \ [base_ptLocal_apex(:,3);1];
                    j3_static_pt_apex4 = base_tx_j3_fixed \ [base_ptLocal_apex(:,4);1];
                    j4_static_pt_apex1 = base_tx_j4_fixed \ [base_ptLocal_apex(:,1);1];
                    
                    % Sovle for joint angles; project corresponding apex
                    % point onto statix XY planes
                    q1 = atan2( j1_static_pt_apex2(2), j1_static_pt_apex2(1) );
                    q2 = atan2( j2_static_pt_apex3(2), j2_static_pt_apex3(1) );
                    q3 = atan2( j3_static_pt_apex4(2), j3_static_pt_apex4(1) );
                    q4 = atan2( j4_static_pt_apex1(2), j4_static_pt_apex1(1) );
                    
                    q = [q1 q2 q3 q4]';
                    
%                     lims=[-200,200;-200,200;-100,300];
%                     items{1}=base_tx_j1_fixed;
%                     items{2}=base_tx_j2_fixed;
%                     items{3}=base_tx_j3_fixed;
%                     items{4}=base_tx_j4_fixed;
%                     items{5}=eye(4,4);
%                     
%                     plotDebug3d(items,10,lims)
%                     
%                     
%                     keyboard
%                     
%                     
%                     
                    
                     
                    % Given apex points in the base frame, solve for kinematic
                    % variables.
                    % joint 1 primary angle correlates to the Y positive apex point
                    % joint 2 primary angle correlates to the X negative apex point
                    % joint 1 complementary angle correlates to the Y negative apex point
                    % joint 2 complementary angle correlates to the X positive apex point
%                     q1 = atan2( base_ptLocal_apex(3,2), base_ptLocal_apex(2,2) );
%                     q2 = atan2( base_ptLocal_apex(3,3), -base_ptLocal_apex(1,3) );
% 
%                     q1c = atan2( base_ptLocal_apex(3,4), -base_ptLocal_apex(2,4) );
%                     q2c = atan2( base_ptLocal_apex(3,1), base_ptLocal_apex(1,1) );
%                     
%                     q(:,ii) = [q1,q2,q1c,q2c]';
                    
                    ikine_soln(ii).q = q;
                    
                    
                    % transform apex points into end effector platform
                    % frame

                    
                    
                    % Solve for angled link rotations; this link rotates
                    % about its Z axis; the X axis is always normal to the
                    % plane defined by the apex points
%                     q1_mid = (obj.alpha/2) - q2c;
%                     q2_mid = (obj.alpha/2) - q1;
%                     q1c_mid = (obj.alpha/2) - q2;
%                     q2c_mid = (obj.alpha/2) - q1c;
%                     
                    
                    link_tx_mid = eye(4,4);
                    link_tx_mid(1:3,1:3) = rotz(pi/2)*rotx(pi/2);
                    link_tx_mid(1:3,4) = [...
                        ((obj.radius/2) / sin(obj.alpha/2)) - ( (obj.ang_len) / tan(obj.alpha/2) );
                        0;
                        -obj.edge/2];
                    
                    mid_tx_top = eye(4,4);
                    mid_tx_top(1:3,1:3) = roty( -obj.alpha ) * rotz( -pi/2 );
                    mid_tx_top(1:3,4) = [...
                        obj.ang_len*cos(obj.alpha)+obj.ang_len,...
                        0,...
                        obj.ang_len*sin(obj.alpha) ]';
            
                    
                    
                    % Create base 4x4 HG TXs for each solution
                    ikine_soln(ii).base_tx_plat = base_txLocal_plat;
                    
                    % Link 1
                    ikine_soln(ii).base_tx_link1_bot = eye(4,4);
                    ikine_soln(ii).base_tx_link1_bot(1:3,1:3) = roty(pi/2)*rotz(pi/2 + q(1) );
                    ikine_soln(ii).base_tx_link1_bot(1:3,4) = [obj.edge/2 0 0]';
                    
                    ikine_soln(ii).base_tx_link1_mid = ikine_soln(ii).base_tx_link1_bot * link_tx_mid;
                    ikine_soln(ii).base_tx_link1_mid(1:3,1:3) = ikine_soln(ii).base_tx_link1_mid(1:3,1:3) * rotz(q3);
%                     
%                     ikine_soln(ii).base_tx_link1_top = ikine_soln(ii).base_tx_link1_mid * mid_tx_top;
%                     ikine_soln(ii).base_tx_link1_top(1:3,1:3) = ikine_soln(ii).base_tx_link1_top(1:3,1:3) *  rotz(dec/2);

           
                    
                    % Link 2
%                     ikine_soln(ii).base_tx_link2_bot = eye(4,4);
%                     ikine_soln(ii).base_tx_link2_bot(1:3,1:3) = rotx(-pi/2)*rotz(pi + q2);
%                     ikine_soln(ii).base_tx_link2_bot(1:3,4) = [0 obj.edge/2 0]';
                    % Link 2 mid
%                     ikine_soln(ii).base_tx_link2_mid = ikine_soln(ii).base_tx_link2_bot * link_tx_mid;
%                     ikine_soln(ii).base_tx_link2_mid(1:3,1:3) = ikine_soln(ii).base_tx_link2_mid(1:3,1:3);% * rotz(q2_mid);
                    % link 2 top
%                     ikine_soln(ii).base_tx_link2_top = ikine_soln(ii).base_tx_link2_mid * mid_tx_top;
%                     ikine_soln(ii).base_tx_link2_top(1:3,1:3) = ikine_soln(ii).base_tx_link2_top(1:3,1:3);% * rotz(q2_mid);
                    
                    
                    
                    % Link 3
%                     ikine_soln(ii).base_tx_link3_bot = eye(4,4);
%                     ikine_soln(ii).base_tx_link3_bot(1:3,1:3) = roty(-pi/2) * rotz(-pi/2 + q1c);
%                     ikine_soln(ii).base_tx_link3_bot(1:3,4) = [-obj.edge/2 0 0]';
                    % link 3 mid
%                     ikine_soln(ii).base_tx_link3_mid = ikine_soln(ii).base_tx_link3_bot * link_tx_mid;
%                     ikine_soln(ii).base_tx_link3_mid(1:3,1:3) = ikine_soln(ii).base_tx_link3_mid(1:3,1:3);% * rotz(q1c_mid);
                    % link 3 top
%                     ikine_soln(ii).base_tx_link3_top = ikine_soln(ii).base_tx_link3_mid * mid_tx_top;
%                     ikine_soln(ii).base_tx_link3_top(1:3,1:3) = ikine_soln(ii).base_tx_link3_top(1:3,1:3);% * rotz(q1c_mid);
                    
                    
                    % Link 4
%                     ikine_soln(ii).base_tx_link4_bot = eye(4,4);
%                     ikine_soln(ii).base_tx_link4_bot(1:3,1:3) = rotx(pi/2) * rotz(q2c);
%                     ikine_soln(ii).base_tx_link4_bot(1:3,4) = [0 -obj.edge/2 0]';
                    % link 4 mid
%                     ikine_soln(ii).base_tx_link4_mid = ikine_soln(ii).base_tx_link4_bot * link_tx_mid;
%                     ikine_soln(ii).base_tx_link4_mid(1:3,1:3) = ikine_soln(ii).base_tx_link4_mid(1:3,1:3);% * rotz(q2c_mid);
                    % link 4 top
%                     ikine_soln(ii).base_tx_link4_top = ikine_soln(ii).base_tx_link4_mid * mid_tx_top;
%                     ikine_soln(ii).base_tx_link4_top(1:3,1:3) = ikine_soln(ii).base_tx_link4_top(1:3,1:3);% * rotz(q2c_mid);
%                     
                    
                    
                    
                    

                    
                end % for
            
            
        end % ikine
        
        function drive()
            
        end % drive
        
        function plot3d( obj, varargin )
            
            % FUNCTION axhan = plot3d( obj, varargin )
            %   Plots a 3d model of the mechansim instance
            %   
            
            soln = [];
            ts = [];
            
            for argn=1:2:nargin-1 %remove the counted obj arg
                
                switch varargin{argn}
                   
                    case 'fighan'
                        obj.fighan = varargin{argn+1};
                        
                    case 'axhan'
                        obj.axhan = varargin{argn+1};
                        
                    case 'soln'
                        soln = varargin{argn+1};
                        
                    case 'ts'
                        ts = varargin{argn+1};
                        
                    otherwise
                        
                    
                end % switch
                    
                
            end % for
            
            %debug
            obj.fighan = [];
            obj.axhan = [];
            
            
            
            if isempty( obj.fighan )
                % Instance figure handle empty; raise a new figure
                obj.fighan = figure( 'name', obj.name,...
                    'numbertitle', 'off', ...
                    'menubar', 'none' );
            end
            
            
            
            if isempty( obj.axhan )
                
                obj.axhan = axes( 'parent', obj.fighan,...
                'projection', 'orthographic', ...
                'units', 'normalized', ....
                'position', [0.1 0.1 0.75 0.75] );
            
            
                obj.axhan.Toolbar.Visible = 'on';
                set( obj.axhan, 'nextplot', 'replacechildren' );
                axis vis3d equal 
                grid on
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                xlim( [-obj.radius,obj.radius].*2 );
                ylim( [-obj.radius,obj.radius].*2);
                zlim( [-obj.radius,obj.radius*3] );
            
                % Since axis was empty, child transforms and child plots
                % have to be drawn
                
                link_len = ((obj.radius/2) / sin(obj.alpha/2)) - ( (obj.ang_len) / tan(obj.alpha/2) );
                
                apex_radius=(obj.radius/2) / sin(obj.alpha/2);
                t=linspace(0,2*pi,1000);
                yy=apex_radius*cos(t);
                zz=apex_radius*sin(t);
                xx=0*t;
                
                world_hgt_base = hgtransform( 'parent', obj.axhan );
                plot3( world_hgt_base, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 3 );
                plot3( world_hgt_base,...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    'linewidth', 1 );
                % plot circle
                plot3( world_hgt_base, ...
                    xx,yy,zz,'k.' );
                    
                
                %debug
                apex_pts = plot3( world_hgt_base, ...
                    soln.base_pt_apex(1,:),soln.base_pt_apex(2,:),soln.base_pt_apex(3,:), 'm*',... 
                    'linewidth', 2 );
                %debug
                
                % Line from base to EE frame; this is a fixed length of
                % obj.radius
                origin_line = plot3( world_hgt_base, ...
                    [0 0], [0 0], [0 0], 'm--', 'linewidth', 1 );
                
                base_hgt_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 3 );
                
                % End effector platform
                base_hgt_plat = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_plat, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 3 );
                plot3( base_hgt_plat,...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    obj.edge/2, 0, 0, 'k*',...
                    -obj.edge/2, 0, 0, 'k*',...
                    0, obj.edge/2, 0, 'k*',...
                    0, -obj.edge/2, 0, 'k*',...
                    'linewidth', 1 );
                % plot circle
                plot3( base_hgt_plat, ...
                    xx,yy,zz,'k.' );
                
                
                % Bottom links
                base_hgt_link1_bot = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link1_bot, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2 );
                plot3( base_hgt_link1_bot, ...
                    [0 link_len link_len], [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                    'linewidth', 1 );
                
                base_hgt_link2_bot = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link2_bot, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link2_bot, ...
%                     [0 link_len link_len], [0 0 0], [0 0 -obj.edge/2], 'k--', ...
%                     'linewidth', 1 );
                
                base_hgt_link3_bot = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link3_bot, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link3_bot, ...
%                     [0 link_len link_len], [0 0 0], [0 0 -obj.edge/2], 'k--', ...
%                     'linewidth', 1 );
                
                base_hgt_link4_bot = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link4_bot, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link4_bot, ...
%                     [0 link_len link_len], [0 0 0], [0 0 -obj.edge/2], 'k--', ...
%                     'linewidth', 1 );
                
                
                base_hgt_link1_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link1_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2 );
                plot3( base_hgt_link1_mid, ...
                    [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                    [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                    'linewidth', 1 );
                
                base_hgt_link2_mid = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link2_mid, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link2_mid, ...
%                     [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
%                     [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
%                     'linewidth', 1 );
                
                base_hgt_link3_mid = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link3_mid, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link3_mid, ...
%                     [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
%                     [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
%                     'linewidth', 1 );
                
                base_hgt_link4_mid = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link4_mid, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link4_mid, ...
%                     [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
%                     [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
%                     'linewidth', 1 );
                
                
                base_hgt_link1_top = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link1_top, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link1_top, ...
%                     [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -link_len], 'k--', ...
%                     'linewidth', 1 );

                
                base_hgt_link2_top = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link2_top, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link2_top, ...
%                     [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -link_len], 'k--', ...
%                     'linewidth', 1 );
%                 
                base_hgt_link3_top = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link3_top, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link3_top, ...
%                     [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -link_len], 'k--', ...
%                     'linewidth', 1 );
                
                base_hgt_link4_top = hgtransform( 'parent', world_hgt_base );
%                 plot3( base_hgt_link4_top, ...
%                     [0 obj.AXLEN], [0 0], [0 0], 'r',...
%                     [0 0], [0 obj.AXLEN], [0 0], 'g',...
%                     [0 0], [0 0], [0 obj.AXLEN], 'b', ...
%                     0,0,0,'k*', 'linewidth', 2 );
%                 plot3( base_hgt_link4_top, ...
%                     [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -link_len], 'k--', ...
%                     'linewidth', 1 );
%                 
                
                
                
                
            end
            
            %debug; add drive controls
            %% UI elements
            pnl = uipanel( 'parent', obj.fighan ,...
                'title', 'Control', ...
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

            
            
            
            for ii=1:length(soln)
                
                axes_update( soln(ii) )
                
             
                % pause
                pause( ts );
                
                drawnow();
                
                
                
            end % for
            
            
            % callback functions
            function resetCallback( varargin )
                
                slider_azimuth.Value = 0;
                slider_declination.Value = 0;

                [~,ss] = obj.ikine(0,0);
                
                axes_update( ss );
            end
            
            function azimuthCallback( varargin )
                azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
                [~,ss] = obj.ikine( slider_declination.Value,slider_azimuth.Value);
                axes_update( ss );
            end
            
            function declinationCallback( varargin )
                declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );
                [~,ss] = obj.ikine( slider_declination.Value,slider_azimuth.Value);
                axes_update( ss );
            end
            
            function axes_update( soln )
                
                apex_pts.XData = soln.base_pt_apex(1,:);
                apex_pts.YData = soln.base_pt_apex(2,:);
                apex_pts.ZData = soln.base_pt_apex(3,:);
                
                
                origin_line.XData = [0 soln.base_tx_plat(1,4)];
                origin_line.YData = [0 soln.base_tx_plat(2,4)];
                origin_line.ZData = [0 soln.base_tx_plat(3,4)];
               
                %update transforms
                base_hgt_plat.Matrix = soln.base_tx_plat;
                base_hgt_mid.Matrix = soln.base_tx_mid;
                base_hgt_link1_bot.Matrix = soln.base_tx_link1_bot;
%                 base_hgt_link2_bot.Matrix = soln.base_tx_link2_bot;
%                 base_hgt_link3_bot.Matrix = soln.base_tx_link3_bot;
%                 base_hgt_link4_bot.Matrix = soln.base_tx_link4_bot;
                
                 base_hgt_link1_mid.Matrix = soln.base_tx_link1_mid;
%                 base_hgt_link2_mid.Matrix = soln.base_tx_link2_mid;
%                 base_hgt_link3_mid.Matrix = soln.base_tx_link3_mid;
%                 base_hgt_link4_mid.Matrix = soln.base_tx_link4_mid;
%                 
%                  base_hgt_link1_top.Matrix = soln.base_tx_link1_top;
%                 base_hgt_link2_top.Matrix = soln.base_tx_link2_top;
%                 base_hgt_link3_top.Matrix = soln.base_tx_link3_top;
%                 base_hgt_link4_top.Matrix = soln.base_tx_link4_top;
                
                
                
            end

            
          
            
        end % plot3d
        
        
    end %methods
    
    
end % classdef