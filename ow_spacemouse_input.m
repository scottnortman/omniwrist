function [dec, az] = ow_spacemouse_input( hDrv, rscale, tscale )
% Assumes a handle to a spacemouse object is passed

% Align with physical mouse on desk
trx = eye(4,4);
trx(1:3,1:3) = rotx(pi/2); % align orientation with mouse on desk

zax = [0 0 1]';


% Y axis is vertical; calc angle from nominal
tx = trx * spacemouse_tx( hDrv, rscale, tscale );

ang = acos( dot( zax, tx(1:3,2) ) / ( norm(zax)*norm(tx(1:3,2)) ) );

if tx(3,2) < 0

    dec = pi/2;
else
    dec = abs( ang );
end

az = atan2( tx(2,2), tx(1,2) );


end