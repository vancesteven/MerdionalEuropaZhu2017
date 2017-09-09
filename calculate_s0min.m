function smin = calculate_s0min(altobeta, fb, fh, smax)
% ARGUMENTS:
% altobeta: beta/alpha ratio;
% fb: geothermal heat flux (W/m2), scalar;
% fh: ice thickness transport (m/s), scalar;
% smax: maximum salinity of the ocean (constrained for MgSO4 and NaCl
%       respectively)
% -------------------------------------------------------------------
% EFFECTS:
% Calculate the minimum average salinity of the ocean by the buoyancy
% requirement.

rho = 1e3;                        % kg/m3, pure water density
rhoi = 920;                       % kg/m3, ice density
cp = 4000;                        % J/(kg k), water heat capacity
smin = altobeta*fb/(rho*cp*rhoi/rho*fh);
if (smin>smax)
    smin = nan;
end
end 