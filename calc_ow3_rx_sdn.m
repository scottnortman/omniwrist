function rx = calc_ow3_rx_sdn( dec, az )

rx = rotx( -dec ) * rotz( -az );


end