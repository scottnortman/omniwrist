function [tx, txMid] = calc_ow3_tx_sdn( dec, az, radius)
% FUNCTION tx = calc_ow3_tx_sdn( dec, az, radius )
%   Given the desired pointing angle of the OmniWrist mechanism as
%   declination and azimuth, returns the corresponding 4x4 homogeneous
%   transform.
%
%   INPUTS
%       
%
%   OUTPUTS
%
%   NOTE
%
%   The simplified kinematic model includes two planar squares, with
%   centers seperated by a constant radius.  The bottom bottom square 
%   represents the fixed base of the mechanism, and
%   the top square represents the moving end effector platform.  Each
%   planar square also includes two coplanar perpendicular lines,
%   intersecting at the center of their respective square.  Said perpendicular lines
%   represent axes of rotation for linkage arms.  
%


    % Calculate the origin point location as a function of 
    % Given the desired declination angle and azimuth, the mechanism end
    % effector platform origin translated along the sphere surface radius, 
    % half of the angle that the platform rotates.  Additionally the Matlab
    % function requires the input as elevation, so we need to convert.
    [x,y,z] = sph2cart( az, (pi/2)-(dec/2), radius );

    % Caclulate the Z axis orientation as a function of azimuth and
    % declination; Note that this angle varies from 0 to pi/2 as the frame
    % origin point subtends an angle of 0 to pi/4, ie half of the
    % orientation angle
    tx=eye(4,4);
    tx(1:3,1:3) = rotz( az ) * roty( dec ) * rotz(-az);
    tx(1:3,4) = [x y z]';
    
    
    
    % calc mid frame
    [xmid,ymid,zmid] = sph2cart( az, (pi/2)-(dec/2), radius/2 );
    
    txMid=eye(4,4);
    txMid(1:3,1:3) = rotz( az ) * roty( dec/2 ) * rotz(-az);
    txMid(1:3,4) = [xmid,ymid,zmid]';
    
    
    
    

end