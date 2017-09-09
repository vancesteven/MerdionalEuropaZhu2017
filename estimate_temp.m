function d_T = estimate_temp(c, u, ks, d, fb)
% ARGUMENTS:
% c: entrainment rate at the base of the freshwater layer;
% u: turbulent velocity (m/s), scalar; 
% ks: effective diffusivity (m2/s), scalar;
% fb: geothermal heat flux (W/m2), scalar.
% --------------------------------------------------------
% EFFECTS:
% Calculate delta T (temperature difference between the deep ocean and the
% layer) based on the heat balance equation of the deep ocean.

cp = 4000;                        % J/(kg k), water heat capacity
rho = 1e3;                        % kg/m3, pure water density

d_T = fb/(rho*cp*(c*u+ks/d));
end

