classdef ow3 < handle
    
    properties
        
        % Instance fixed parameters
        radius = [];        % Distance between base and platform centers
        edge = [];          % Edge length of base and platform squares
        alpha = [];         % radians, total angle between upper/lower arms
        ang_len =[];        % Length of edge of angle link offset, optional
        link_len=[];        % Length of linkage from joint axis to apex
        name = 'ow3';       % Instance name
        
        
        fighan = [];
        axhan = [];
       
        
        AXLEN = [];
        

        % Homogeneous transforms of the actual instance state
        
        % to do : add assertions in code
        declination_min = 0;
        declination_max = pi/2;
        azimuth_min =  -pi/2;
        azimuth_max = pi/2;
               
        
        
        base_declination_plat = 0; %add these to code
        base_azimuth_plat = 0;
        

        
        base_hgt_plat = eye(4,4);
        base_hgt_cen = eye(4,4);
        
        base_hgt_j1_fixed = eye(4,4);
        base_hgt_j2_fixed = eye(4,4);
        base_hgt_j3_fixed = eye(4,4);
        base_hgt_j4_fixed = eye(4,4);
        
        base_hgt_link1_bot = eye(4,4);
        base_hgt_link2_bot = eye(4,4);
        base_hgt_link3_bot = eye(4,4);
        base_hgt_link4_bot = eye(4,4);
        
        base_hgt_link1_mid = eye(4,4);
        base_hgt_link2_mid = eye(4,4);
        base_hgt_link3_mid = eye(4,4);
        base_hgt_link4_mid = eye(4,4);
        
        base_hgt_link1_top = eye(4,4);
        base_hgt_link2_top = eye(4,4);
        base_hgt_link3_top = eye(4,4);
        base_hgt_link4_top = eye(4,4);
        
        link_tx_ang = eye(4,4);
        mid_tx_top = eye(4,4);
        
        
        base_ptx_ptA = eye(4,4);
        base_ptx_ptB = eye(4,4);
        base_ptx_ptC = eye(4,4);
        base_ptx_ptD = eye(4,4);
          
        
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
            
            
            
            % Geometric link length from point of rotation to apex point
            obj.link_len = ((obj.radius/2) / sin(obj.alpha/2));
            
            % Define fixed frames 
            obj.link_tx_ang = eye(4,4);
            obj.link_tx_ang(1:3,1:3) = rotz(pi/2)*rotx(pi/2);
            obj.link_tx_ang(1:3,4) = [...
                obj.link_len - ( (obj.ang_len) / tan(obj.alpha/2) );
                0;
                -obj.edge/2];

            obj.mid_tx_top = eye(4,4);
            obj.mid_tx_top(1:3,1:3) = roty( -obj.alpha ) * rotz( -pi/2 );
            obj.mid_tx_top(1:3,4) = [...
                obj.ang_len*cos(obj.alpha)+obj.ang_len,...
                0,...
                obj.ang_len*sin(obj.alpha) ]';
            
            % Fixed joint coordinate frames, j1

            obj.base_hgt_j1_fixed(1:3,1:3) = roty(pi/2) * rotz(pi/2);
            obj.base_hgt_j1_fixed(1:3,4) = [obj.edge/2, 0, 0]';
            % Fixed joint coordinate frames, j2

            obj.base_hgt_j2_fixed(1:3,1:3) = rotx(-pi/2) * rotz(pi);
            obj.base_hgt_j2_fixed(1:3,4) = [0 obj.edge/2, 0]';
            % Fixed joint coordinate frames, j3

            obj.base_hgt_j3_fixed(1:3,1:3) = roty(-pi/2) * rotz(-pi/2);
            obj.base_hgt_j3_fixed(1:3,4) = [-obj.edge/2, 0, 0]';
            % Fixed joint coordinate frames, j4

            obj.base_hgt_j4_fixed(1:3,1:3) = rotx(pi/2);
            obj.base_hgt_j4_fixed(1:3,4) = [0, -obj.edge/2, 0]';
            
            
            obj.drive( 0, 0 );
            
            
        end %ow3
        
        %{
        
        function ji = jin( obj, declination, azimuth )
            % function ji = jin( obj, declination, azimuth )
            %   Numerically calculates an approximation of the inverse 
            %   analytical jacobian about the 
            %   specified task space operating point.      
            %
            %   INPUTS
            %       declination,    scalar, [0, pi/2]
            %       azimuth,        scalar, [-pi/2, pi/2]
            %
            %   OUTPUTS
            %       ji,             2x6 inverse jacobian matrix
            %
            %   NOTE
            %       The inverse jacobian ji is a 2x6 matrix to map the 6x1
            %       spatial twist to joint space velocity:
            %
            %           qd = ji * xd
            %
            %       
            
            % This is done
            % numerically by linear approximation about a given
            % operation x.
            %
            %       qd = ji * xd
            %
            %   Given the individual inverse jacobian matrix elements
            %   dqn / dx  
            %
            %   [   ji(1,1) ji(1,2) ji(1,3) ji(1,4) ji(1,5) ji(1,6);
            %       ji(2,1) ji(2,2) ji(2,3) ji(2,4) ji(2,5) ji(2,6) ];
            %
            %   Rearrange into an A * x = b form by setting 
            %
            %       x = ji'
            %
            %   And each row of the regression matrix A is 
            %
            %       A = [dx dy dz drx dry drz]
            % 
            %   and each row of the b matrix is
            %
            %       b = [dq1 dq2]
            %   
            %
            %  So we can numerically solve for x with
            %
            %    [dx dy dz drx dry drz] * x = [dq1 dq2]
            %
            %   Where x = ji'
            %
            %   Differential task space values are calculated by
            %   calculating the homogeneous transform differences
            % tdx = inv(tx0) * tx1
            
   
            
        end % 
        
        %}
        
        function [declination, azimuth, base_tx_plat] = fkine( obj, q )
            % function [declination, azimuth, base_tx_plat] = fkine( obj, q )
            %   Calculates the forward kinematic solution of the mechanism
            %   given the specified joint values q.
            %
            %   INPUTS
            %       q,  2 x 1 x n, column vector of joint values
            %
            %   OUTPUTS
            
            %
            %   NOTE
            %       Solution based on the paper:
            %
            %       "Geometric Approach for Kinematic Analysis of a Class 
            %       of 2-DOF Rotational Parallel Manipulators" by Dong Xin,
            %       Chen Bing, and Zong Guanghua, Chinese Journal of 
            %       Mechanical Engineering, Oct 2011
            %
            
            % Each soln is 2x1 column vector
            [~,nc] = size( q );
            
            declination = zeros( nc,1 );
            azimuth = zeros( nc, 1 );
            base_tx_plat = zeros( 4,4,nc );     
            
            for ii=1:nc
            
                % Calculate apex points in local frames given q
                j1_pt_apex = rotz( q(1,ii) ) * [obj.link_len 0 -obj.edge/2]';
                j2_pt_apex = rotz( q(2,ii) ) * [obj.link_len 0 -obj.edge/2]';

                %Get points as referenced in the paper
                base_pt_J1 = obj.base_hgt_j1_fixed * [j1_pt_apex; 1];
                base_pt_J1 = base_pt_J1(1:3,1);
                base_pt_J2 = obj.base_hgt_j2_fixed * [j2_pt_apex; 1];
                base_pt_J2 = base_pt_J2(1:3,1);

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
                declination(ii) = pi - 2*atan( base_pt_Op_fkine(3,1) / sqrt( base_pt_Op_fkine(1,1)^2 + base_pt_Op_fkine(2,1)^2 ) );

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

                azimuth(ii) = atan2( yv, xv );

                [x,y,z] = sph2cart( azimuth, (pi/2)-(declination/2), obj.radius );
                base_tx_plat(:,:,ii) = eye(4,4);
                base_tx_plat(1:3,1:3,ii) = rotz( azimuth ) * roty( declination ) * rotz( -azimuth );
                base_tx_plat(1:3,4,ii) = [x y z]';
            
            end % for
            
        end % fkine
        
        %{
        function jf = jfn( obj, declination, azimuth )
            % function jf = jfn( obj, declination, azimuth )
            %   Calculates the numerical forward jacobian at the given pose
            %   using a numerical finite difference approach
            
            %   The forward jacobian maps the joint rates with respect to
            %   time, qd, to the task space rates xd with respect to time:
            %
            %       xd = jfn * qd
            %
            %   In the case of the OmniwristIII mechanism, since there are
            %   only two controlled joints q1 and q1, the shape of the
            %   matrices are
            %
            %       xd  => 6 x 1
            %       [   dx;
            %           dy;
            %           dz;
            %           dwx;
            %           dwy;
            %           dwz; ];
            %           
            %
            %       jfn => 6 x 2
            %       [   dx/dq1  dx/dq2;
            %           dy/dq1  dy/dq2;
            %           dz/dq1  dz/dq2;
            %           dwx/dq1 dwx/dq2;
            %           dwy/dq1 dwy/dq2;
            %           dwz/dq1 dwz/dq2; ];
            %       
            %       qd  => 2 x 1
            %       [   dq1;
            %           dq2; ];
            %
            
            
            
            
            
        end % jfn
        %}
        
        
        function [q, base_tx_plat, base_tx_cen, base_pt_apex] = ikine2( obj, declination, azimuth )
            
            ld = length( declination );
            la = length( azimuth );
            
            q = [];
            base_tx_plat = [];

            if ld ~= la
                warning('Declination and azimuth vector lengths must match' );
                return;
            end
            
            
            
            q = zeros( 4, ld );
            base_pt_apex = zeros(3,4,ld);
    
            base_tx_plat=zeros(4,4,ld);  
            
            % todo move these fixed transforms to object properties

            
            % Solve for each pair
            for ii=1:ld
                
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
                % orientation angle, along the sphere surface per above;
                % add in the parasitic motion about the Z axis as a
                % function of azimuth
                base_tx_plat(:,:,ii) = eye(4,4);
                base_tx_plat(1:3,1:3,ii) = rotz( az ) * roty( dec ) * rotz( -az );
                base_tx_plat(1:3,4,ii) = [x y z]';
                
                % Calculate mid frame
                [xm, ym, zm] = sph2cart( az, (pi/2)-(dec/2), obj.radius/2 );
                base_tx_cen = eye(4,4);
                base_tx_cen(1:3,1:3) = rotz( az ) * roty( dec/2 ) * rotz( -az );
                base_tx_cen(1:3,4) = [xm, ym, zm]';
                
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
                %   solve the two remaining unknown parameters.
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
                thetaX_q1q3 = cir_zero_comp_solve( obj.link_len, ...    % Circle radius
                    base_tx_plat(1,4,ii), ...                              % X component of position
                    base_tx_plat(1,3,ii), ...                              % X component of platform Z axis
                    base_tx_plat(1,2,ii) ) ;                               % X component of platform Y axis


                if any( isnan( thetaX_q1q3 ) ) || any( isinf(thetaX_q1q3) )
                   % If result returns NaN, assume circles are co-planar
                   % and therefore solve with the 2D solution

                   % Solution is in the YZ plane                        
                   cen1 = [0 0]';
                   cen2 = [base_tx_plat(3,4,ii),base_tx_plat(2,4,ii)]';
                   r1 = obj.link_len;
                   r2 = obj.link_len;
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

                    q13xyz_A =     base_tx_plat(1:3,4,ii) ...
                            +   obj.link_len*base_tx_plat(1:3,3,ii)*cos( theta_q1 ) ... % platform Z axis is 0 ref for theta
                            +   obj.link_len*base_tx_plat(1:3,2,ii)*sin( theta_q1 );
                    q13xyz_B =     base_tx_plat(1:3,4,ii) ...
                            +   obj.link_len*base_tx_plat(1:3,3,ii)*cos( theta_q3 ) ... % platform Z axis is 0 ref for theta
                            +   obj.link_len*base_tx_plat(1:3,2,ii)*sin( theta_q3 );

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
                thetaY_q2q4 = cir_zero_comp_solve( obj.link_len, ... % Circle radius
                    base_tx_plat(2,4,ii), ...                      % Y component of position
                    base_tx_plat(2,1,ii), ...                      % Y component of platform X axis
                    base_tx_plat(2,3,ii) );                        % Y component of platform Z axis

                if any( isnan( thetaY_q2q4 ) ) || any( isinf(thetaY_q2q4 ) )
                    % If result returns NaN, assume circles are co-planar
                    % and therefore solve with the 2D solution

                    % Solution is in the XZ plane                        
                    cen1 = [0 0]';
                    cen2 = [base_tx_plat(1,4,ii),base_tx_plat(3,4,ii)]';
                    r1 = obj.link_len;
                    r2 = obj.link_len;
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

                    q24xyz_A =     base_tx_plat(1:3,4,ii) ...
                            +   obj.link_len*base_tx_plat(1:3,1,ii)*cos( theta_q2 ) ...
                            +   obj.link_len*base_tx_plat(1:3,3,ii)*sin( theta_q2 );
                    q24xyz_B =     base_tx_plat(1:3,4,ii) ...
                            +   obj.link_len*base_tx_plat(1:3,1,ii)*cos( theta_q4 ) ... 
                            +   obj.link_len*base_tx_plat(1:3,3,ii)*sin( theta_q4 );

                    % Q2 is neg x, q4 is pos x
                    if q24xyz_A(1) > q24xyz_B(1)
                        q4xyz = q24xyz_A;
                        q2xyz = q24xyz_B;
                    else
                        q4xyz = q24xyz_B;
                        q2xyz = q24xyz_A;
                    end

                end %q2q4
               
                % Given the apex points, solve for each link joint
                % angle; transform into local frames
                j1_pt_q1xyz = obj.base_hgt_j1_fixed \ [q1xyz;1];
                j2_pt_q2xyz = obj.base_hgt_j2_fixed \ [q2xyz;1];
                j3_pt_q3xyz = obj.base_hgt_j3_fixed \ [q3xyz;1];
                j4_pt_q4xyz = obj.base_hgt_j4_fixed \ [q4xyz;1];

                % Sovle for joint angles; project corresponding apex
                % point onto statix XY planes
                q1 =  atan2( j1_pt_q1xyz(2), j1_pt_q1xyz(1) );
                q2 =  atan2( j2_pt_q2xyz(2), j2_pt_q2xyz(1) );
                q3 =  atan2( j3_pt_q3xyz(2), j3_pt_q3xyz(1) );
                q4 =  atan2( j4_pt_q4xyz(2), j4_pt_q4xyz(1) );
                
                % NoteL  q1,q2 are actuated joints; q3,q4 are passive
                q(:,ii) = [q1 q2 q3 q4]';
                
                base_pt_apex(:,1,ii) = q1xyz;
                base_pt_apex(:,2,ii) = q2xyz;
                base_pt_apex(:,3,ii) = q3xyz;
                base_pt_apex(:,4,ii) = q4xyz;
                
                
            end % for
            
        end % ikine2
  
        function drive( obj, declination, azimuth )
            % function drive( obj, declination, azimuth )
            % Updates the moving frames of the object instance
            
            
            ld = length( declination );
            la = length( azimuth );
            
            
            assert( ld == la );
            
            for ii=1:ld 
                
                dec = declination( ii );
                az = azimuth( ii );
                
%                 assert( ( obj.declination_min <= dec) && ( dec <= obj.declination_max ) );
%                 assert( ( obj.azimuth_min <= az ) && (az <= obj.azimuth_max) );
                
                obj.base_declination_plat = dec;
                obj.base_azimuth_plat = az;
        
            
                % values
                %[q, base_tx_plat, base_tx_cen, apex] = ikine2( obj, declination, azimuth )
                [q, obj.base_hgt_plat, obj.base_hgt_cen, base_pt_apex] = obj.ikine2( declination, azimuth );
                
                obj.base_ptx_ptA(1:3,4) = base_pt_apex(:,1);
                obj.base_ptx_ptB(1:3,4) = base_pt_apex(:,2);
                obj.base_ptx_ptC(1:3,4) = base_pt_apex(:,3);
                obj.base_ptx_ptD(1:3,4) = base_pt_apex(:,4);

                % Link 1 moving frames
                obj.base_hgt_link1_bot = obj.base_hgt_j1_fixed;
                obj.base_hgt_link1_bot(1:3,1:3) = obj.base_hgt_link1_bot(1:3,1:3) * rotz( q(1) );
                % Define nominal mid frame
                obj.base_hgt_link1_mid = obj.base_hgt_link1_bot * obj.link_tx_ang;
                % Transform platform frame origin into this frame
                link1_mid_pt_plat = obj.base_hgt_link1_mid \ obj.base_hgt_plat(:,4);
                % Calculate mid, top rotation (they are equal)
                qmt1 = pi/2 - atan2( link1_mid_pt_plat(1), link1_mid_pt_plat(2) );
                obj.base_hgt_link1_mid(1:3,1:3) = obj.base_hgt_link1_mid(1:3,1:3) * rotz( qmt1 );
                % Top frame
                obj.base_hgt_link1_top = obj.base_hgt_link1_mid * obj.mid_tx_top;
                obj.base_hgt_link1_top(1:3,1:3) = obj.base_hgt_link1_top(1:3,1:3) * rotz( qmt1 );


                % Link 2 moving frames
                obj.base_hgt_link2_bot = obj.base_hgt_j2_fixed;
                obj.base_hgt_link2_bot(1:3,1:3) = obj.base_hgt_link2_bot(1:3,1:3) * rotz( q(2) );
                % Define nominal mid frame
                obj.base_hgt_link2_mid = obj.base_hgt_link2_bot * obj.link_tx_ang;
                % Transform platform frame origin into this frame
                link2_mid_pt_plat = obj.base_hgt_link2_mid \ obj.base_hgt_plat(:,4);
                % Calculate mid, top rotation (they are equal)
                qmt2 = pi/2 - atan2( link2_mid_pt_plat(1), link2_mid_pt_plat(2) );
                obj.base_hgt_link2_mid(1:3,1:3) = obj.base_hgt_link2_mid(1:3,1:3) * rotz( qmt2 );
                % Top frame
                obj.base_hgt_link2_top = obj.base_hgt_link2_mid * obj.mid_tx_top;
                obj.base_hgt_link2_top(1:3,1:3) = obj.base_hgt_link2_top(1:3,1:3) * rotz( qmt2 );

                % Link 3 moving frames
                obj.base_hgt_link3_bot = obj.base_hgt_j3_fixed;
                obj.base_hgt_link3_bot(1:3,1:3) = obj.base_hgt_link3_bot(1:3,1:3) * rotz( q(3) );
                % Define nominal mid frame
                obj.base_hgt_link3_mid = obj.base_hgt_link3_bot * obj.link_tx_ang;
                % Transform platform frame origin into this frame
                link3_mid_pt_plat = obj.base_hgt_link3_mid \ obj.base_hgt_plat(:,4);
                % Calculate mid, top rotation (they are equal)
                qmt3 = pi/2 - atan2( link3_mid_pt_plat(1), link3_mid_pt_plat(2) );
                obj.base_hgt_link3_mid(1:3,1:3) = obj.base_hgt_link3_mid(1:3,1:3) * rotz( qmt3 );
                % Top frame
                obj.base_hgt_link3_top = obj.base_hgt_link3_mid * obj.mid_tx_top;
                obj.base_hgt_link3_top(1:3,1:3) = obj.base_hgt_link3_top(1:3,1:3) * rotz( qmt3 );

                % Link 4 moving frames
                obj.base_hgt_link4_bot = obj.base_hgt_j4_fixed;
                obj.base_hgt_link4_bot(1:3,1:3) = obj.base_hgt_link4_bot(1:3,1:3) * rotz( q(4) );
                % Define nominal mid frame
                obj.base_hgt_link4_mid = obj.base_hgt_link4_bot * obj.link_tx_ang;
                % Transform platform frame origin into this frame
                link4_mid_pt_plat = obj.base_hgt_link4_mid \ obj.base_hgt_plat(:,4);
                % Calculate mid, top rotation (they are equal)
                qmt4 = pi/2 - atan2( link4_mid_pt_plat(1), link4_mid_pt_plat(2) );
                obj.base_hgt_link4_mid(1:3,1:3) = obj.base_hgt_link4_mid(1:3,1:3) * rotz( qmt4 );
                % Top frame
                obj.base_hgt_link4_top = obj.base_hgt_link4_mid * obj.mid_tx_top;
                obj.base_hgt_link4_top(1:3,1:3) = obj.base_hgt_link4_top(1:3,1:3) * rotz( qmt4 );
                
                % If the instance axis handle is not empty, assume we have
                % to update the homogeneous transform matricies in the axis
                try 
                    if isvalid( obj.axhan )

                        % Find object HGT properties
                        pp = properties( obj );
                        pl = pp(contains( pp ,{'hgt', 'ptx' } )  );

                        % Update HGT matricies
                        for jj=1:length( pl )

                            hgt_tag = pl{jj};
                            mtx = eval( ['obj.', hgt_tag] );
                            set( findobj( obj.axhan, 'tag', hgt_tag ), 'Matrix', mtx ) ;

                        end                      
                        
                    end
                    
                catch ME
                    disp(ME);
                    keyboard
                end
                               
            end % for
            
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
                xlim( [-obj.radius,obj.radius].*2 );
                ylim( [-obj.radius,obj.radius].*2);
                zlim( [-obj.radius,obj.radius*2] );
            
                % Since axis was empty, child transforms and child plots
                % have to be drawn
                
                hgtransform( 'parent', obj.axhan, 'tag', 'world_hgt_base' );
                plot3( findobj( obj.axhan, 'tag', 'world_hgt_base' ), ...
                    [0 obj.AXLEN], [0 0], [0 0], 'r',...
                    [0 0], [0 obj.AXLEN], [0 0], 'g',...
                    [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                    0,0,0,'k*', 'linewidth', 3, 'tag', 'frame' );
                % Generate fixed base lines
                plot3( findobj( obj.axhan, 'tag', 'world_hgt_base' ),...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    'linewidth', 1, 'tag', 'line' );
                % Generate points for a circle
                apex_radius=(obj.radius/2) / sin(obj.alpha/2);
                t=linspace(0,2*pi,1000);
                yy=apex_radius*cos(t);
                zz=apex_radius*sin(t);
                xx=0*t;
                % plot circles
                plot3( findobj( obj.axhan, 'tag', 'world_hgt_base' ), ...
                    xx,yy,zz,'k.', 'tag', 'circle' );
                plot3( findobj( obj.axhan, 'tag', 'world_hgt_base' ), ...
                    yy,xx,zz,'k.', 'tag', 'circle' );
                
                
                
                
                
                % Create needed coordinate frames to plot; this corresponds
                % to object properties that include the '_hgt_' substring
                pp = properties( obj );
                pl = pp(contains( pp ,'hgt')  );
                for jj=1:length( pl )

                    try
                        
                        hgt_tag = pl{jj};
                        hgtransform( 'parent', findobj( obj.axhan, 'tag', 'world_hgt_base' ), 'tag', hgt_tag );
                        
                        
                        mtx = eval( ['obj.', hgt_tag] );
                        set( findobj( obj.axhan, 'tag', hgt_tag ), 'Matrix', mtx ) ;
                        
                        % Plot coordinate frame
                        plot3( findobj( obj.axhan, 'tag', hgt_tag ), ...
                            [0 obj.AXLEN], [0 0], [0 0], 'r',...
                            [0 0], [0 obj.AXLEN], [0 0], 'g',...
                            [0 0], [0 0], [0 obj.AXLEN], 'b', ...
                            0,0,0,'k*', 'linewidth', 3, 'tag', 'frame' );

                        % Plot bottom linkage lines
                        if contains( hgt_tag, '_bot' )
                            plot3( findobj( obj.axhan, 'tag', hgt_tag ), ...
                                [0 (obj.link_len - (obj.ang_len/tan(obj.alpha/2))) (obj.link_len - (obj.ang_len/tan(obj.alpha/2)))], ...
                                [0 0 0], [0 0 -obj.edge/2], 'k--', ...
                                'linewidth', 1, 'tag', 'line' );
                        end
                        
                        % Plot middle angle linkage lines
                        if contains( hgt_tag, '_mid' )
                        
                            plot3( findobj( obj.axhan, 'tag', hgt_tag ),...
                                [0 obj.ang_len obj.ang_len*cos(obj.alpha)+obj.ang_len], ...
                                [0 0 0], [0 0 obj.ang_len*sin(obj.alpha)], 'k--', ...
                                'linewidth', 1, 'tag', 'line' );
                        end
                        
                        % Plot top linkage lines
                        if contains( hgt_tag, '_top' )
                            plot3( findobj( obj.axhan, 'tag', hgt_tag ), ...
                                [0 obj.edge/2 obj.edge/2], [0 0 0], ...
                                [0 0 -(obj.link_len - (obj.ang_len/tan(obj.alpha/2))) ], 'k--', ...
                                'linewidth', 1, 'tag', 'line' ); 
                        end
                        
                        if contains( hgt_tag, '_cen' )
                            plot3( findobj( obj.axhan, 'tag', hgt_tag ), ...
                                [0 0], [0 0], [-obj.radius/2 obj.radius/2], ...
                                'm--', 'linewidth', 1, 'tag', 'line' );
                        end
                    
                    
                    catch ME
                       disp(ME);
                        keyboard 
                    end
                    
                end % for
                

                % Generate fixed base lines
                plot3( findobj( obj.axhan, 'tag', 'base_hgt_plat' ),...
                    [0 obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 -obj.edge/2], [0 0], [0 0], 'k--',...
                    [0 0], [0 obj.edge/2], [0 0], 'k--',...
                    [0 0], [0 -obj.edge/2], [0 0], 'k--',...
                    'linewidth', 1, 'tag', 'line' );
                plot3( findobj( obj.axhan, 'tag', 'base_hgt_plat' ),...
                    obj.edge/2, 0, 0, 'k*',...
                    -obj.edge/2, 0, 0, 'k*',...
                    0, obj.edge/2,0, 'k*',...
                    0, -obj.edge/2, 0, 'k*',...
                    'linewidth', 3, 'tag', 'point' );
                % plot circles
                plot3( findobj( obj.axhan, 'tag', 'base_hgt_plat' ), ...
                    xx,yy,zz,'k.', 'tag', 'circle' );
                plot3( findobj( obj.axhan, 'tag', 'base_hgt_plat' ), ...
                    yy,xx,zz,'k.', 'tag', 'circle' );
                
                % Circle intersection points
                pp = properties( obj );
                pl = pp(contains( pp ,'ptx')  );
                pls = {'r*', 'g*', 'b*', 'c*'};
                for jj=1:length( pl )
                    ptx_tag = pl{jj};
                    hgtransform( 'parent', findobj( obj.axhan, 'tag', 'world_hgt_base' ), 'tag', ptx_tag );
                    plot3( findobj( obj.axhan, 'tag', ptx_tag ),...
                        0,0,0, pls{jj},  'linewidth', 2, 'tag', 'point' );
                    mtx = eval( ['obj.', ptx_tag] );
                    set( findobj( obj.axhan, 'tag', ptx_tag ), 'Matrix', mtx ) ;
                end
                            
                
                
                
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
                'value', obj.base_azimuth_plat, ...
                'position', [0 0.025 1, 0.4],...
                'callback', @azimuthCallback );

            %Declination slider
            slider_declination = uicontrol( sld_pnl,...
                'style', 'slider', ...
                'units', 'normalized',...
                'min', 0, 'max', pi/2,... %min was 0
                'value', obj.base_declination_plat, ...
                'position', [0 0.5 1, 0.4],...
                'callback', @declinationCallback );

            % Azimuth text value
            azimuth_text = uicontrol( sld_pnl, ...
                'style', 'text', ...
                'enable', 'inactive', ...
                'string', sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi ), ...
                'units', 'normalized', ...
                'fontweight', 'bold',...
                'position', [0 0.425, 1, 0.05]);


            % Declination text value
            declination_text = uicontrol( sld_pnl, ...
                'style', 'text', ...
                'enable', 'inactive', ...
                'string', sprintf( 'Azi: [%0.2f]°', slider_declination.Value*180/pi ), ...
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
  
                slider_azimuth.Value = 0;
                slider_declination.Value = 0;
                
                azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
                declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );

                axes_update();

            end
            
            function azimuthCallback( varargin )
                
                azimuth_text.String = sprintf( 'Azi: [%0.2f]°', slider_azimuth.Value*180/pi );
                
                if abs(slider_declination.Value) <= 100*eps
                   slider_declination.Value = 0;
                end
                
                if abs(slider_azimuth.Value) <= 100*eps
                   slider_azimuth.Value = 0;
                end
                axes_update();
            end
            
            function declinationCallback( varargin )
                
                
                declination_text.String = sprintf( 'Dec: [%0.2f]°', slider_declination.Value*180/pi );
                
                if abs(slider_declination.Value) <= 100*eps
                   slider_declination.Value = 0;
                end
                
                if abs(slider_azimuth.Value) <= 100*eps
                   slider_azimuth.Value = 0;
                end
                
                axes_update();
            end
            
            function axes_update( )
                
                % how to update text values when calling 'drive'
                % externally?
                
                obj.drive( slider_declination.Value, slider_azimuth.Value );
                
            end

            
          
            
        end % plot3d
        
        
    end %methods
    
    
end % classdef