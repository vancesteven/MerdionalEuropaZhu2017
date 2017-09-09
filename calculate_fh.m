function fh = calculate_fh(Ts_e, Ts_p, dFocn)
% ARGUMENTS:
% Ts_e: the surface temperature at equator;
% Ts_p: the surface temperature at pole;
% dFocn: delta Focn, the equator-to-pole difference in ocean-ice heat flux
% ------------------------------------------------------------------------
% EEFECTS:
% Calcuate ice thickness transport (fh, m/s) based on the thickness balance
% equation.

ki = 2;                                  % thermal conductivity of ice, W/m/K
h0 = 1e4;                                % equilibrium ice thickness, m
Li = 3.3e8;                              % J/m3, latent heat of fusion
delta_Ts = Ts_e - Ts_p; 
fh = ki*delta_Ts/(2*h0*Li) + dFocn/2/Li;
end