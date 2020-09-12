classdef ow3 < handle
    
    properties
        
        % Instance constants
        radius = [];        % Distance between kinematic centers
        edge = [];          % Edge length of base and platform squares
        alpha = [];         % radians, total angle between upper/lower arms
        ang_len =[];        % Length of edge of angle link (optional)
        link_rad=[];
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
        
        tx_list = { };
        
        
        
        
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
            
            
            obj.link_rad = ((obj.radius/2) / sin(obj.alpha/2));
            
            
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
                
                                
                len = length(declination);                
                q = zeros(4, len );
                ikine_soln = struct( ...
                    'q', cell(1,len),...
                    'q1xyz', cell(1,len),...
                    'q2xyz', cell(1,len),...
                    'q3xyz', cell(1,len),...
                    'q4xyz', cell(1,len),...
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
                    base_tx_plat=eye(4,4);
                    base_tx_plat(1:3,1:3) = rotz( az ) * roty( dec ) * rotz( -az );
                    base_tx_plat(1:3,4) = [x y z]';
                    
                    % Inverse kinematic solution requires the calculation
                    % of the intersection points of two circles in 3d
                    % space.  We can define an arbitrary circle
                    % parametrically with the following:
                    %
                    %   P = C + r*cos(theta)*U + r*sin(theta)*V
                    %
                    %   where
                    %
                    %   P, 3x1, point on circle
                    %   C, 3x1, center point of circle   
                    %   r, scalar, circle radius
                    %   U, 3x1, unit vector perpendicular to circle axis
                    %   V, 3x1, unit vector, V = cross( axis, U )
                    %   theta, scalar parameter, [0,2*pi]
                    
                    % For each arm, we can simplify the solution given that
                    % the intersection points will always lie in either the
                    % YZ or XZ plane of the base coordinate system
                    %   q1 => x = 0, solution in YZ plane, y => positive
                    %   q2 => y = 0, solution in XZ plane, x => negative
                    %   q3 => x = 0, solution in YZ plane, y => negative
                    %   q4 => y = 0, solution in XZ plane, x => positive
                    %
                    %   Since we know the above components, we can solve
                    %   for the parametric value theta, and then calculate
                    %   the remaining unknown components using the
                    %   parametric equation with the solved theta value to
                    %   calculate the two remaining unknown parameters.
                    % 
                    %   Note this assumption is not valid when said circles
                    %   are coplanar; in this case, an alternative 2DOF
                    %   solution is used (only needed when declination is
                    %   0).
                    
                    % cir_zero_comp_solve(r,c0,u0,v0)
                    % q1, q3 x=0 in the base plane, with the x component 0
                    % This solves for the parametric theta value of a
                    % circle in the platform YZ plane , centered at the
                    % platform origin, intersecting with the
                    % base frame YZ plane (ie, the theta value when the x
                    % component is 0)
                    thetaX_q1q3 = cir_zero_comp_solve( obj.link_rad, ...    % Circle radius
                        base_tx_plat(1,4), ...                              % X component of position
                        base_tx_plat(1,3), ...                              % X component of platform Z axis
                        base_tx_plat(1,2) ) ;                               % X component of platform Y axis
                                          

                    if any( isnan( thetaX_q1q3 ) ) || any( isinf(thetaX_q1q3) )
                       % If result returns NaN, assume circles are co-planar
                       % and therefore solve with the 2D solution
                        
                       % Solution is in the YZ plane                        
                       cen1 = [0 0]';
                       cen2 = [base_tx_plat(3,4),base_tx_plat(2,4)]';
                       r1 = obj.link_rad;
                       r2 = obj.link_rad;
                       % Solution is [z,y]
                       [p1,p2] = cir_cir_int( cen1, r1, cen2, r2 );
                       
                       % for q1 pick solution with positive y coordinate
                       % and q3 has a 
                       if p1(2) > 0
                           q1p = p1;
                           q3p = p2;
                       else
                           q1p = p2;
                           q3p = p1;
                       end
                           
                       % Both solutions have X = 0
                       q1xyz = [0 q1p(2) q1p(1)]';
                       q3xyz = [0 q3p(2) q3p(1)]';
                      
                    else
                        
                        theta_q1 = real( thetaX_q1q3(1) );
                        theta_q3 = real( thetaX_q1q3(2) );

                        q13xyz_A =     base_tx_plat(1:3,4) ...
                                +   obj.link_rad*base_tx_plat(1:3,3)*cos( theta_q1 ) ... % platform Z axis is 0 ref for theta
                                +   obj.link_rad*base_tx_plat(1:3,2)*sin( theta_q1 );
                        q13xyz_B =     base_tx_plat(1:3,4) ...
                                +   obj.link_rad*base_tx_plat(1:3,3)*cos( theta_q3 ) ... % platform Z axis is 0 ref for theta
                                +   obj.link_rad*base_tx_plat(1:3,2)*sin( theta_q3 );
                            
                        % q1 has posY, q3 has neg Y
                        if q13xyz_A(2) > q13xyz_B(2)
                            q1xyz = q13xyz_A;
                            q3xyz = q13xyz_B;
                        else
                            q1xyz = q13xyz_B;
                            q3xyz = q13xyz_A;
                        end
                        
                    end % q1q3
                    
                    % q2, q4 y=0 in the base plane, with the y component 0
                    % This solves for the parametric theta value of a
                    % circle in the platform XZ plane, centered at the
                    % platform origin, intersecting with the base frame XZ
                    % plane (ie, the theta value when the y component is 0)
                    thetaY_q2q4 = cir_zero_comp_solve( obj.link_rad, ... % Circle radius
                        base_tx_plat(2,4), ...                      % Y component of position
                        base_tx_plat(2,1), ...                      % Y component of platform X axis
                        base_tx_plat(2,3) );                        % Y component of platform Z axis

                    if any( isnan( thetaY_q2q4 ) ) || any( isinf(thetaY_q2q4 ) )
                        % If result returns NaN, assume circles are co-planar
                        % and therefore solve with the 2D solution

                        % Solution is in the XZ plane                        
                        cen1 = [0 0]';
                        cen2 = [base_tx_plat(1,4),base_tx_plat(3,4)]';
                        r1 = obj.link_rad;
                        r2 = obj.link_rad;
                        % Solution is [x, z]
                        [p1,p2] = cir_cir_int( cen1, r1, cen2, r2 );

                        %  Determine which point is for q2 or q4 based on X
                        %  value
                        if p1(1) > 0
                            % p1 is for q4, p2 is q2
                           q2p = p2;
                           q4p = p1;
                        else
                           q2p = p1;
                           q4p = p2;
                        end

                        % Both solutions have y = 0
                        q2xyz = [q2p(1) 0 q2p(2)]';
                        q4xyz = [q4p(1) 0 q4p(2)]';

                    else

                        theta_q2 = real( thetaY_q2q4(1) );
                        theta_q4 = real( thetaY_q2q4(2) );
                        
                        q24xyz_A =     base_tx_plat(1:3,4) ...
                                +   obj.link_rad*base_tx_plat(1:3,1)*cos( theta_q2 ) ...
                                +   obj.link_rad*base_tx_plat(1:3,3)*sin( theta_q2 );
                        q24xyz_B =     base_tx_plat(1:3,4) ...
                                +   obj.link_rad*base_tx_plat(1:3,1)*cos( theta_q4 ) ... 
                                +   obj.link_rad*base_tx_plat(1:3,3)*sin( theta_q4 );
                        
                        % Q2 is neg x, q4 is pos x
                        if q24xyz_A(1) > q24xyz_B(1)
                            q4xyz = q24xyz_A;
                            q2xyz = q24xyz_B;
                        else
                            q4xyz = q24xyz_B;
                            q2xyz = q24xyz_A;
                        end
                    
                    end %q2q4
                   
                    
                    % Create static base coordinate frames for each joint
                    % note q1 => +X, 
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
                    
                    
                    % Given the apex points, solve for each link joint
                    % angle; transform into local frames
                    base_pt_q1xyz = base_tx_j1_fixed \ [q1xyz;1];
                    base_pt_q2xyz = base_tx_j2_fixed \ [q2xyz;1];
                    base_pt_q3xyz = base_tx_j3_fixed \ [q3xyz;1];
                    base_pt_q4xyz = base_tx_j4_fixed \ [q4xyz;1];
                    
                    % Sovle for joint angles; project corresponding apex
                    % point onto statix XY planes
                    q1 =  atan2( base_pt_q1xyz(2), base_pt_q1xyz(1) );
                    q2 =  atan2( base_pt_q2xyz(2), base_pt_q2xyz(1) );
                    q3 =  atan2( base_pt_q3xyz(2), base_pt_q3xyz(1) );
                    q4 =  atan2( base_pt_q4xyz(2), base_pt_q4xyz(1) );
                    
                    % Define fixed frames 
                    link_tx_mid = eye(4,4);
                    link_tx_mid(1:3,1:3) = rotz(pi/2)*rotx(pi/2);
                    link_tx_mid(1:3,4) = [...
                        obj.link_rad - ( (obj.ang_len) / tan(obj.alpha/2) );
                        0;
                        -obj.edge/2];
                    
                    mid_tx_top = eye(4,4);
                    mid_tx_top(1:3,1:3) = roty( -obj.alpha ) * rotz( -pi/2 );
                    mid_tx_top(1:3,4) = [...
                        obj.ang_len*cos(obj.alpha)+obj.ang_len,...
                        0,...
                        obj.ang_len*sin(obj.alpha) ]';
  
                    % Create base 4x4 HG TXs for each solution
                    ikine_soln(ii).base_tx_plat = base_tx_plat;
                    
                    % Link 1
                    ikine_soln(ii).base_tx_link1_bot = base_tx_j1_fixed;
                    ikine_soln(ii).base_tx_link1_bot(1:3,1:3) = ikine_soln(ii).base_tx_link1_bot(1:3,1:3) * rotz( q1 );
                    % Define nominal mid frame
                    ikine_soln(ii).base_tx_link1_mid = ikine_soln(ii).base_tx_link1_bot * link_tx_mid;
                    % Transform platform frame origin into this frame
                    link1_mid_pt_plat = ikine_soln(ii).base_tx_link1_mid \ base_tx_plat(:,4);
                    % Calculate mid, top rotation (they are equal)
                    qmt1 = pi/2 - atan2( link1_mid_pt_plat(1), link1_mid_pt_plat(2) );
                    ikine_soln(ii).base_tx_link1_mid(1:3,1:3) = ikine_soln(ii).base_tx_link1_mid(1:3,1:3) * rotz( qmt1 );
                    % Top frame
                    ikine_soln(ii).base_tx_link1_top = ikine_soln(ii).base_tx_link1_mid * mid_tx_top;
                    ikine_soln(ii).base_tx_link1_top(1:3,1:3) = ikine_soln(ii).base_tx_link1_top(1:3,1:3) * rotz( qmt1 );
                    
                    
                    % Link 2
                    ikine_soln(ii).base_tx_link2_bot = base_tx_j2_fixed;
                    ikine_soln(ii).base_tx_link2_bot(1:3,1:3) = ikine_soln(ii).base_tx_link2_bot(1:3,1:3) * rotz( q2 );
                    % Define nominal mid frame
                    ikine_soln(ii).base_tx_link2_mid = ikine_soln(ii).base_tx_link2_bot * link_tx_mid;
                    % Transform platform frame origin into this frame
                    link2_mid_pt_plat = ikine_soln(ii).base_tx_link2_mid \ base_tx_plat(:,4);
                    % Calculate mid, top rotation (they are equal)
                    qmt2 = pi/2 - atan2( link2_mid_pt_plat(1), link2_mid_pt_plat(2) );
                    ikine_soln(ii).base_tx_link2_mid(1:3,1:3) = ikine_soln(ii).base_tx_link2_mid(1:3,1:3) * rotz( qmt2 );
                    % Top frame
                    ikine_soln(ii).base_tx_link2_top = ikine_soln(ii).base_tx_link2_mid * mid_tx_top;
                    ikine_soln(ii).base_tx_link2_top(1:3,1:3) = ikine_soln(ii).base_tx_link2_top(1:3,1:3) * rotz( qmt2 );
                    
                    % Link 3
                    ikine_soln(ii).base_tx_link3_bot = base_tx_j3_fixed;
                    ikine_soln(ii).base_tx_link3_bot(1:3,1:3) = ikine_soln(ii).base_tx_link3_bot(1:3,1:3) * rotz( q3 );
                    % Define nominal mid frame
                    ikine_soln(ii).base_tx_link3_mid = ikine_soln(ii).base_tx_link3_bot * link_tx_mid;
                    % Transform platform frame origin into this frame
                    link3_mid_pt_plat = ikine_soln(ii).base_tx_link3_mid \ base_tx_plat(:,4);
                    % Calculate mid, top rotation (they are equal)
                    qmt3 = pi/2 - atan2( link3_mid_pt_plat(1), link3_mid_pt_plat(2) );
                    ikine_soln(ii).base_tx_link3_mid(1:3,1:3) = ikine_soln(ii).base_tx_link3_mid(1:3,1:3) * rotz( qmt3 );
                    % Top frame
                    ikine_soln(ii).base_tx_link3_top = ikine_soln(ii).base_tx_link3_mid * mid_tx_top;
                    ikine_soln(ii).base_tx_link3_top(1:3,1:3) = ikine_soln(ii).base_tx_link3_top(1:3,1:3) * rotz( qmt3 );
                    
                    % Link 4
                    ikine_soln(ii).base_tx_link4_bot = base_tx_j4_fixed;
                    ikine_soln(ii).base_tx_link4_bot(1:3,1:3) = ikine_soln(ii).base_tx_link4_bot(1:3,1:3) * rotz( q4 );
                    % Define nominal mid frame
                    ikine_soln(ii).base_tx_link4_mid = ikine_soln(ii).base_tx_link4_bot * link_tx_mid;
                    % Transform platform frame origin into this frame
                    link4_mid_pt_plat = ikine_soln(ii).base_tx_link4_mid \ base_tx_plat(:,4);
                    % Calculate mid, top rotation (they are equal)
                    qmt4 = pi/2 - atan2( link4_mid_pt_plat(1), link4_mid_pt_plat(2) );
                    ikine_soln(ii).base_tx_link4_mid(1:3,1:3) = ikine_soln(ii).base_tx_link4_mid(1:3,1:3) * rotz( qmt4 );
                    % Top frame
                    ikine_soln(ii).base_tx_link4_top = ikine_soln(ii).base_tx_link4_mid * mid_tx_top;
                    ikine_soln(ii).base_tx_link4_top(1:3,1:3) = ikine_soln(ii).base_tx_link4_top(1:3,1:3) * rotz( qmt4 );                  
         
                    
                    ikine_soln(ii).q = [q1 q2 q3 q4]';
                    
                    ikine_soln(ii).q1xyz = q1xyz;
                    ikine_soln(ii).q2xyz = q2xyz;
                    ikine_soln(ii).q3xyz = q3xyz;
                    ikine_soln(ii).q4xyz = q4xyz;
                    
                    ikine_soln(ii).base_tx_plat = base_tx_plat;
                    
                    q = ikine_soln(ii).q;

                end % for
            
            
        end % ikine
        
        function fkine( obj, q )
            
            
            
        end % fkine
        
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
                'projection', 'perspective', ... %orthographic or perspective
                'units', 'normalized', ....
                'position', [0.1 0.1 0.75 0.75] );
            
            
                obj.axhan.Toolbar.Visible = 'on';
                set( obj.axhan, 'nextplot', 'replacechildren' );
                axis vis3d equal 
                grid on
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                xlim( [-obj.radius,obj.radius].*1.5 );
                ylim( [-obj.radius,obj.radius].*1.5);
                zlim( [-obj.radius,obj.radius*2] );
            
                % Since axis was empty, child transforms and child plots
                % have to be drawn
                
                % Generate points for a circle
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
                    0,0,0,'k*', 'linewidth', 3, 'tag', 'frame' );
                plot3( world_hgt_base,...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    'linewidth', 1, 'tag', 'line' );
%                 plot circles
                plot3( world_hgt_base, ...
                    xx,yy,zz,'k.', 'tag', 'circle' );
                plot3( world_hgt_base, ...
                    yy,xx,zz,'k.', 'tag', 'circle' );
                    
                
                % Line from base to EE frame; this is a fixed length of
                % obj.radius
                origin_line = plot3( world_hgt_base, ...
                    [0 0], [0 0], [0 0], 'm--', 'linewidth', 1, 'tag', 'line' );
                
                q1_pt = plot3( world_hgt_base,...
                    soln.q1xyz(1), soln.q1xyz(2), soln.q1xyz(3), 'r*', ...
                    'linewidth', 2, 'tag', 'point' );
                q2_pt = plot3( world_hgt_base,...
                    soln.q2xyz(1), soln.q2xyz(2), soln.q2xyz(3), 'g*', ...
                    'linewidth', 2, 'tag', 'point' );
                q3_pt = plot3( world_hgt_base,...
                    soln.q3xyz(1), soln.q3xyz(2), soln.q3xyz(3), 'b*', ...
                    'linewidth', 2, 'tag', 'point'  );
                q4_pt = plot3( world_hgt_base,...
                    soln.q4xyz(1), soln.q4xyz(2), soln.q4xyz(3), 'c*', ...
                    'linewidth', 2, 'tag', 'point'  );
                
                

                
                % Forward kinematic debug; based on the paper:
                % "Geometric Approach for Kinematic Analysis of a Class of
                % 2-DOF Rotational Parallel Manipulators" by Dong Xin, Chen
                % Bing, and Zong Guanghua, Chinese Journal of Mechanical
                % Engineering, Oct 2011
                % Draw line fom apex points
                qall = [soln.q1xyz,soln.q2xyz,soln.q3xyz,soln.q4xyz,soln.q1xyz];
                qall_pt = plot3( world_hgt_base,...
                    qall(1,:), qall(2,:),qall(3,:), 'k', ...
                    'linewidth', 2, 'tag', 'circle' );
                
                %Get points as referenced in the paper
                base_pt_J1 = soln.q1xyz;
                base_pt_J2 = soln.q2xyz;
                base_pt_O=[0 0 0]';
                base_pt_J12 = base_pt_J2-base_pt_J1;
                J12_len = norm( base_pt_J12 );
                base_pt_J12_unit = base_pt_J12 / J12_len;
                base_pt_M = base_pt_J1 + base_pt_J12_unit*J12_len/2;
                % Create a coordinate frame with origin = M, Z axis along
                % MJ2, and X along MO
                base_pt_origin = base_pt_M;
                base_pt_zaxis = base_pt_J2-base_pt_M;
                base_pt_zaxis = base_pt_zaxis / norm(base_pt_zaxis);
                base_pt_xaxis = base_pt_O - base_pt_M;
                base_pt_xaxis = base_pt_xaxis / norm( base_pt_xaxis );
                base_pt_yaxis = cross( base_pt_zaxis, base_pt_xaxis );
                base_pt_yaxis = base_pt_yaxis / norm( base_pt_yaxis );
                base_tx_midJ1J2 = eye(4,4);
                base_tx_midJ1J2(1:3,1) = base_pt_xaxis;
                base_tx_midJ1J2(1:3,2) = base_pt_yaxis;
                base_tx_midJ1J2(1:3,3) = base_pt_zaxis;
                base_tx_midJ1J2(1:3,4) = base_pt_origin;
                
                % Calculate the base_pt_O in the mid frame
                midJ1J2_pt_O = base_tx_midJ1J2 \ [base_pt_O;1];
                
                % Calculate the vector length from the mid point frame
                % origin to the pt_O; used to solve the forward kinematics
                len_ptO = norm(midJ1J2_pt_O(1:3,1));
                
                % Need to calculate Op in the midJ1J2 frame; since we know
                % points M (origin of midJ1J2), O and Op are in a plane,
                % and we know the fixed length of O to Op, we can
                % calculate the needed rotation angle rz  to rotate
                % vector MO to MOp (from law of cosines)
                rz = acos( (len_ptO^2 + len_ptO^2 - obj.radius^2) / (2*len_ptO^2) );
                
                % Rotate vector MO by rz to get MOp
                midJ1J2_pt_Op = rotz( rz ) * midJ1J2_pt_O(1:3,1);
                
                % Transform back into base frame to get final solution
                base_pt_Op_fkine = base_tx_midJ1J2 * [midJ1J2_pt_Op;1];
                
                % Given vector base_OOp, calculate corresponding declination and
                % azimuth values
                dec = pi - 2*atan( base_pt_Op_fkine(3,1) / sqrt( base_pt_Op_fkine(1,1)^2 + base_pt_Op_fkine(2,1)^2 ) );
                
                % Handle special case close to 0
                if abs( base_pt_Op_fkine( 2,1 )) < 100*eps
                    yv = 0;
                else
                    yv = base_pt_Op_fkine( 2,1 );
                end
                
                if abs( base_pt_Op_fkine( 1,1 ) ) < 100*eps
                    
                    xv = base_pt_Op_fkine( 1,1 );
                else
                    xv = base_pt_Op_fkine( 1,1 );
                end
                    
                az = atan2( yv, xv );
                
                
                
                
               
                
                
  
                
                base_hgt_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                base_hgt_mid.Matrix = base_tx_midJ1J2;
                
                
                
                
                
                
                
                % End effector platform
                base_hgt_plat = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_plat, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 3, 'tag', 'frame' );
                plot3( base_hgt_plat,...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    obj.edge/2, 0, 0, 'k*',...
                    -obj.edge/2, 0, 0, 'k*',...
                    0, obj.edge/2, 0, 'k*',...
                    0, -obj.edge/2, 0, 'k*',...
                    'linewidth', 1, 'tag', 'line' );
                % plot circles
                plot3( base_hgt_plat, ...
                    xx,yy,zz,'k.', 'tag', 'circle' );
                plot3( base_hgt_plat, ...
                    yy,xx,zz,'k.', 'tag', 'circle' );
                
                

                % Bottom links
                base_hgt_link1_bot = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link1_bot, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link1_bot, ...
                    [0 (obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) (obj.link_rad - (obj.ang_len/tan(obj.alpha/2)))], ...
                    [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link2_bot = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link2_bot, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link2_bot, ...
                    [0 (obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) (obj.link_rad - (obj.ang_len/tan(obj.alpha/2)))], ...
                    [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link3_bot = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link3_bot, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link3_bot, ...
                    [0 (obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) (obj.link_rad - (obj.ang_len/tan(obj.alpha/2)))], ...
                    [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link4_bot = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link4_bot, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link4_bot, ...
                    [0 (obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) (obj.link_rad - (obj.ang_len/tan(obj.alpha/2)))], ...
                    [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                
                base_hgt_link1_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link1_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link1_mid, ...
                    [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                    [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link2_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link2_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link2_mid, ...
                    [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                    [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link3_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link3_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame'  );
                plot3( base_hgt_link3_mid, ...
                    [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                    [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link4_mid = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link4_mid, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link4_mid, ...
                    [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                    [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                
                base_hgt_link1_top = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link1_top, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link1_top, ...
                    [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -(obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) ], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );

                
                base_hgt_link2_top = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link2_top, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link2_top, ...
                    [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -(obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) ], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
%                 
                base_hgt_link3_top = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link3_top, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link3_top, ...
                    [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -(obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) ], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                base_hgt_link4_top = hgtransform( 'parent', world_hgt_base );
                plot3( base_hgt_link4_top, ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 2, 'tag', 'frame' );
                plot3( base_hgt_link4_top, ...
                    [0 obj.edge/2 obj.edge/2], [0 0 0], [0 0 -(obj.link_rad - (obj.ang_len/tan(obj.alpha/2))) ], 'k--', ...
                    'linewidth', 1, 'tag', 'line' );
                
                
                
                
                
            end
            
            %debug; add drive controls
            %% UI elements
            sld_pnl = uipanel( 'parent', obj.fighan ,...
                'title', 'Slide', ...
                'units', 'normalized',...
                'position', [0.775 0 0.2 0.75] );

            reset_btn = uicontrol( sld_pnl, ...
                'style', 'pushbutton', ...
                'units', 'normalized', ...
                'position', [0 0.95, 1, 0.05], ...
                'string', 'RESET', ...
                'callback', @resetCallback );

            %Azimuth slider
            slider_azimuth = uicontrol( sld_pnl,...
                'style', 'slider', ...
                'units', 'normalized',...
                'min', -pi, 'max', pi,...
                'value', 0, ...
                'position', [0 0.025 1, 0.4],...
                'callback', @azimuthCallback );

            %Declination slider
            slider_declination = uicontrol( sld_pnl,...
                'style', 'slider', ...
                'units', 'normalized',...
                'min', 0, 'max', pi/2,... %min was 0
                'value', 0, ...
                'position', [0 0.5 1, 0.4],...
                'callback', @declinationCallback );

            % Azimuth text value
            azimuth_text = uicontrol( sld_pnl, ...
                'style', 'text', ...
                'enable', 'inactive', ...
                'string', 'Azimuth' , ...
                'units', 'normalized', ...
                'fontweight', 'bold',...
                'position', [0 0.425, 1, 0.05]);


            % Declination text value
            declination_text = uicontrol( sld_pnl, ...
                'style', 'text', ...
                'enable', 'inactive', ...
                'string', 'Declination' , ...
                'units', 'normalized', ...
                'fontweight', 'bold',...
                'position', [0 0.9, 1, 0.05] );
            
            ctl_pnl = uipanel( 'parent', obj.fighan ,...
                'title', 'Control', ...
                'units', 'normalized',...
                'position', [0.775 0.775 0.2 0.2] );
            
            % visual element check boxes
            cir_ctl = uicontrol( ctl_pnl, ...
                'style', 'checkbox', ...
                'units', 'normalized', ...
                'string', 'Circles', ...
                'position', [0 0.85 0.2 0.2],...
                'value', 1, ...
                'callback', @plotFeatureCallback );
            line_ctl = uicontrol( ctl_pnl, ...
                'style', 'checkbox', ...
                'units', 'normalized', ...
                'string', 'Lines', ...
                'position', [0 0.65 0.2 0.2],...
                'value', 1, ...
                'callback', @plotFeatureCallback );
            frm_ctl = uicontrol( ctl_pnl, ...
                'style', 'checkbox', ...
                'units', 'normalized', ...
                'string', 'Frames', ...
                'position', [0 0.45 0.2 0.2],...
                'value', 1, ...
                'callback', @plotFeatureCallback );
            pts_ctl = uicontrol( ctl_pnl, ...
                'style', 'checkbox', ...
                'units', 'normalized', ...
                'string', 'Points', ...
                'position', [0 0.25 0.2 0.2],...
                'value', 1, ...
                'callback', @plotFeatureCallback );
            stl_ctl = uicontrol( ctl_pnl, ...
                'style', 'checkbox', ...
                'units', 'normalized', ...
                'string', '3D Models', ...
                'position', [0 0.05 0.2 0.2],...
                'value', 1, ...
                'callback', @plotFeatureCallback );
            
            

            
            
            
            for ii=1:length(soln)
                
                axes_update( soln(ii) )
                
             
                % pause
                pause( ts );
                
                drawnow();
                
                
                
            end % for
            
            
            % callback functions
            
            function plotFeatureCallback( varargin )
                
                %circles
                if cir_ctl.Value
                   set(findobj( obj.axhan, 'tag', 'circle' ), 'visible', 'on' )
                else
                    set(findobj( obj.axhan, 'tag', 'circle' ), 'visible', 'off' )
                end
                
                %line
                if line_ctl.Value
                   set(findobj( obj.axhan, 'tag', 'line' ), 'visible', 'on' )
                else
                    set(findobj( obj.axhan, 'tag', 'line' ), 'visible', 'off' )
                end
                
                %frames
                if frm_ctl.Value
                   set(findobj( obj.axhan, 'tag', 'frame' ), 'visible', 'on' )
                else
                    set(findobj( obj.axhan, 'tag', 'frame' ), 'visible', 'off' )
                end
                
                %points
                if pts_ctl.Value
                   set(findobj( obj.axhan, 'tag', 'point' ), 'visible', 'on' )
                else
                    set(findobj( obj.axhan, 'tag', 'point' ), 'visible', 'off' )
                end
                
                %stls
                if stl_ctl.Value
                   set(findobj( obj.axhan, 'tag', 'stl' ), 'visible', 'on' )
                else
                    set(findobj( obj.axhan, 'tag', 'stl' ), 'visible', 'off' )
                end
                
            end
            
            
            function resetCallback( varargin )
                
                clc();
                
                slider_azimuth.Value = 0;
                slider_declination.Value = 0;
                
                azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
                declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );

                [~,ss] = obj.ikine(0,0);
                
                axes_update( ss );
            end
            
            function azimuthCallback( varargin )
                
                clc();
                
                azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
                
                if abs(slider_declination.Value) <= 100*eps
                   slider_declination.Value = 0;
                end
                
                if abs(slider_azimuth.Value) <= 100*eps
                   slider_azimuth.Value = 0;
                end
                
                
                [~,ss] = obj.ikine( slider_declination.Value,slider_azimuth.Value );
                axes_update( ss );
            end
            
            function declinationCallback( varargin )
                
                clc();
                
                declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );
                
                if abs(slider_declination.Value) <= 100*eps
                   slider_declination.Value = 0;
                end
                
                if abs(slider_azimuth.Value) <= 100*eps
                   slider_azimuth.Value = 0;
                end
                
                
                [~,ss] = obj.ikine( slider_declination.Value,slider_azimuth.Value);
                
                
                axes_update( ss );
            end
            
            function axes_update( soln )
                
                % debug
                soln.q

                
                origin_line.XData = [0 soln.base_tx_plat(1,4)];
                origin_line.YData = [0 soln.base_tx_plat(2,4)];
                origin_line.ZData = [0 soln.base_tx_plat(3,4)];
                
                q1_pt.XData = soln.q1xyz(1);
                q1_pt.YData = soln.q1xyz(2);
                q1_pt.ZData = soln.q1xyz(3);
                
                q2_pt.XData = soln.q2xyz(1);
                q2_pt.YData = soln.q2xyz(2);
                q2_pt.ZData = soln.q2xyz(3);
                
                q3_pt.XData = soln.q3xyz(1);
                q3_pt.YData = soln.q3xyz(2);
                q3_pt.ZData = soln.q3xyz(3);
                
                q4_pt.XData = soln.q4xyz(1);
                q4_pt.YData = soln.q4xyz(2);
                q4_pt.ZData = soln.q4xyz(3);
                
                qall = [soln.q1xyz,soln.q2xyz,soln.q3xyz,soln.q4xyz,soln.q1xyz];
                
                qall_pt.XData = qall(1,:);
                qall_pt.YData = qall(2,:);
                qall_pt.ZData = qall(3,:);
                
                
                %%%fkine debug
    
                %Get points as referenced in the paper
                base_pt_J1 = soln.q1xyz;
                base_pt_J2 = soln.q2xyz;
                base_pt_O=[0 0 0]';
                base_pt_J12 = base_pt_J2-base_pt_J1;
                J12_len = norm( base_pt_J12 );
                base_pt_J12_unit = base_pt_J12 / J12_len;
                base_pt_M = base_pt_J1 + base_pt_J12_unit*(J12_len/2);
                % Create a coordinate frame with origin = M, Z axis along
                % MJ2, and X along MO
                base_pt_origin = base_pt_M;
                base_pt_zaxis = base_pt_J2-base_pt_M;
                base_pt_zaxis = base_pt_zaxis / norm(base_pt_zaxis);
                base_pt_xaxis = base_pt_O - base_pt_M;
                base_pt_xaxis = base_pt_xaxis / norm( base_pt_xaxis );
                base_pt_yaxis = cross( base_pt_zaxis, base_pt_xaxis );
                base_pt_yaxis = base_pt_yaxis / norm( base_pt_yaxis );
                base_tx_midJ1J2 = eye(4,4);
                base_tx_midJ1J2(1:3,1) = base_pt_xaxis;
                base_tx_midJ1J2(1:3,2) = base_pt_yaxis;
                base_tx_midJ1J2(1:3,3) = base_pt_zaxis;
                base_tx_midJ1J2(1:3,4) = base_pt_origin;
                
                midJ1J2_pt_O = base_tx_midJ1J2 \ [base_pt_O;1];
                len_ptO = norm(midJ1J2_pt_O(1:3,1));
                
                % Need to calculate Op in the midJ1J2 frame; since we know
                % points M (origin of midJ1J2), 0 and Op are in a plane,
                % and we know the fixed distance of O to Op, we can
                % calculate the needed rotation angle alpha to rotate
                % vector MO to MOp
                rz = acos( (len_ptO^2 + len_ptO^2 - obj.radius^2) / (2*len_ptO^2) );
                
                midJ1J2_pt_Op = rotz( rz ) * midJ1J2_pt_O(1:3,1);
                base_pt_Op_fkine = base_tx_midJ1J2 * [midJ1J2_pt_Op;1];
                
                dec = pi - 2*atan( base_pt_Op_fkine(3,1) / sqrt( base_pt_Op_fkine(1,1)^2 + base_pt_Op_fkine(2,1)^2 ) );
                % Handle special case close to 0
                if abs( base_pt_Op_fkine( 2,1 )) < 100*eps
                    yv = 0;
                else
                    yv = base_pt_Op_fkine( 2,1 );
                end
                
                if abs( base_pt_Op_fkine( 1,1 ) ) < 100*eps
                    
                    xv = base_pt_Op_fkine( 1,1 );
                else
                    xv = base_pt_Op_fkine( 1,1 );
                end
                    
                az = atan2( yv, xv );
                
                
                fkine = [...
                    soln.base_tx_plat(1:3,4),...
                    base_pt_Op_fkine(1:3,1)]
                
                cmd =  [...
                    slider_declination.Value, dec;...
                    slider_azimuth.Value, az ] 
               
                
                
                
                %%%
                

                
               
                %update transforms
                base_hgt_plat.Matrix = soln.base_tx_plat;
               base_hgt_link1_bot.Matrix = soln.base_tx_link1_bot;
                base_hgt_link2_bot.Matrix = soln.base_tx_link2_bot;
                base_hgt_link3_bot.Matrix = soln.base_tx_link3_bot;
                base_hgt_link4_bot.Matrix = soln.base_tx_link4_bot;
                
                 base_hgt_link1_mid.Matrix = soln.base_tx_link1_mid;
                base_hgt_link2_mid.Matrix = soln.base_tx_link2_mid;
                base_hgt_link3_mid.Matrix = soln.base_tx_link3_mid;
                base_hgt_link4_mid.Matrix = soln.base_tx_link4_mid;
                
                 base_hgt_link1_top.Matrix = soln.base_tx_link1_top;
                base_hgt_link2_top.Matrix = soln.base_tx_link2_top;
                base_hgt_link3_top.Matrix = soln.base_tx_link3_top;
                base_hgt_link4_top.Matrix = soln.base_tx_link4_top;
                
                
                
            end

            
          
            
        end % plot3d
        
        
    end %methods
    
    
end % classdef