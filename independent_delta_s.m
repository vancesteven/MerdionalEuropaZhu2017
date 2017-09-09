function [delta_s, ri, c] =independent_delta_s(u, d, fh, s0, beta, ks)
% ARGUMENTS:
% u: turbulent velocity (m/s), scalar; 
% d: depth of freshwater layer (m), scalar;
% fh: ice thickness transport (m/s), scalar;
% s0: average salinity of the ocean (psu), scalar;
% beta: haline contraction coefficient (psu^-1);
% ks: effective diffusivity (m2/s), scalar;
% -----------------------------------------------------------------------
% RETURNS:
% delta_s: the two positive roots of delta S (salinity difference between 
%          the deep ocean (S0) and the freshwater layer);
% ri: Richardson number at delta_s;
% c: entrainment rate at delta_s;
% -----------------------------------------------------------------------
% EFFECTS:
% Calculate salinity contrast between the deep ocean and the freshwater
% based on the salinity balance equation of the layer (analytically
% converted to a third-order polynomial equation of delta S^0.5 with no
% assumptions or estimations).
% Return only the positive roots of delta S; the negative root is disgarded. 

g = 1.3;                          % m/s2, gravitational acceleration on Europa
rho = 1e3;                        % kg/m3, pure water density
rhoi = 920;                       % kg/m3, ice density

p(1) = (rhoi/rho)*fh + ks/d;
p(2) = 0;
p(3) = -s0*(rhoi/rho)*fh;
p(4) = 1.5*u^4*(g*beta*d)^(-3/2);
temp = roots(p);
delta_s = (temp(2:3)).^2;
ri = (g*beta*delta_s*d)/u^2;       
c = 1.5*ri.^(-3/2);
end