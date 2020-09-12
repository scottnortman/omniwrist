function rx = calc_ow3_rx( dec, az )
% FUNCTION x = calc_ow3_rx( dec, az )
%   Calculates the 3x3 rotation matrix of the end effector platform given
%   the declination and azimuth values.
%
%   Taken from 
%   https://www.researchgate.net/profile/Jozef_Sofka/publication/
%   224381503_Integrated_approach_to_electromechanical_design_of_a_
%   digitally_controlled_high_precision_actuator_for_aerospace_applications
%   /links/554ca93d0cf21ed2135be40a.pdf

%   https://www.researchgate.net/profile/Jozef_Sofka/publication/224381503_Integrated_approach_to_electromechanical_design_of_a_digitally_controlled_high_precision_actuator_for_aerospace_applications/links/554ca93d0cf21ed2135be40a.pdf

    cd = cos( dec );
    ca = cos( az );
    sd = sin( dec );
    sa = sin( az );
    
    rx = [      cd      0       -sd;    ...
                sd*sa   -ca     -cd*sa; ...
                sd*ca   -sa     -cd*ca; ];

end