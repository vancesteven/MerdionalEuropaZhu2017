function [s2, s2_atmin, d_T, dT, dmin_tur_num, dmax_num] = calculate_dT_at_dmin_num(u, s0, d, beta, altobeta,fh, ks, fb)
% ARGUMENTS:
% u: turbulent velocity (m/s), could be an array; 
% s0: average salinity of the ocean (psu), should be a scalar;
% d: depth of freshwater layer (m), could be an array;
% beta: haline contraction coefficient (psu^-1);
% altobeta: beta/alpha ratio;
% fh: ice thickness transport (m/s), scalar;
% ks: effective diffusivity (m2/s), scalar;
% fb: geothermal heat flux (W/m2), scalar;
%----------------------------------------------------------------
% RETURNS:
% s2: delta S (salinity difference between the deep ocean and the
%     freshwater layer), s2(u, d);
% s2_atmin: delta S at the minimum depth (dmin);
% d_T: delta T (temperature difference between the deep ocean and the
%      layer)at each (u, d), d_T(u, d);
% dT: delta T at dmin, dT(dmin, u);
% dmin_tur_num: dmin (minimum depth of the freshwater layer by both the real
%                     root and turbulent regime requirements),
%                     dmin_tur_num[1,length(u)];
% dmax_num: dmax (minimum depth of the freshwater layer by the suppressing
%                 effect requirement);
%-----------------------------------------------------------------
% EFFECTS:
% Based on the salinity balance of the freshwater layer and the heat
% balance of the deep ocean, calculate delta T, delta S and the critical
% depths of the freshwater layer all numerically (without using analytical
% solutions).

delta_s = zeros(2, size(d, 2), size(u, 2));      % hold two roots of the salinity equation
ri = zeros(2, size(d, 2), size(u, 2));           % Richardson number
c1 = zeros(2, size(d, 2), size(u, 2));           % entrainment rate at the base of the freshwater layer
dmin = zeros(1, length(u));                      % minimum depth by the real root requirement

% Calculate delta S and dmin for real root criterion numerically
for j=1:size(u,2)
    firstreal = 1;
    for i=1:size(d,2)
        [delta_s0, ri0, c10] = independent_delta_s(u(j), d(i), fh, s0, beta, ks);
        if (isreal(delta_s0))
            if (firstreal)
                dmin(j) = d(i);
                firstreal = 0;
            end 
            delta_s(:,i,j) = delta_s0;
            ri(:,i,j) = ri0;
            c1(:,i,j) = c10;
        end
    end
end
dmin(dmin==0)=nan;                                   % no dmin exists if no real root for delta S for all u 
s2 = squeeze(permute(delta_s(2,:,:),[1,3,2]));       % extract the larger root of delta S, which is in a turbulent regime
s2(s2==0) = NaN;                                     % no real root for delta S

% Estimate the temperature of the deep ocean (assuming T of the layer is Tf)
c1 = squeeze(permute(c1(2,:,:),[1,3,2]));            % c1(u, d)
if (length(u)==1)
    c1 = permute(c1(:,:),[2,1]);
    s2 = permute(s2(:,:),[2,1]);
end
c1(c1==0) = NaN; 

% Find dmin and dmax of the layer
dmin_tur_num = zeros(1, length(u));                  % miniumum depth by the turbulent regime requirement
dmax_num = zeros(1, length(u));                      % maximum depth by the suppressing effect requirement
ind_dmin = zeros(1, length(u));                      % index of d that gives dmin_tur_num 
for i=1:length(u)
    % find numerical dmin of the turbulent regime
    for j = 1:length(d)
        if (c1(i,j)*u(i) > ks/d(j))                  % condition of the turbulent regime
            dmin_tur_num(i) = d(j);
            ind_dmin(i) = j;                         % record the index of d that gives dmin_tur to retrieve associated dT
            break;
        end        
    end
    
    % find numerical dmax of the suppressing effect
    for j = 1:length(d)
        if (c1(i,j)>1e-3)                            % condion of the suppressing effect is c << 1e-3
            if (j>1)
                dmax_num(i) = d(j-1);
            else
                dmax_num(i) = nan;
            end
            break;
        else
            if (j==length(d))
                dmax_num(i) = 1e5;                  % at this u, if all d's satisfy the requirement, dmax = total ocean depth
            end
        end
    end
end
dmin_tur_num(dmin_tur_num == 0) = nan;
dmax_num(dmax_num == 0) = nan;

% Calculate temperature contrast at each u and d
d_T = zeros(size(u, 2), size(d, 2));
for i=1:size(u,2)
    for j=1:size(d,2)
        d_T(i,j) = estimate_temp(c1(i,j), u(i), ks, d(j),fb);
    end 
end

% Calculate alpha*delta T/(beta*delta S)
ratio = altobeta*d_T./s2;
ratio(ratio>1) = NaN;

% Determine the candidate for the final minimum depth (the larger one of dmin and dmin_tur)
for i=1:size(u, 2)
    if (i<=size(dmin, 2))
        if (dmin_tur_num(i)< dmin(i))               
            dmin_tur_num(i)= dmin(i);
        end
    else
        dmin_tur_num(i)=NaN;
    end    

    % Find the index of the first d that satisfies the buoyant requirement (ratio<1)
    j = 1;                                          
    while isnan(ratio(i, j))
        j = j+1;
        if j>length(d)
            break;
        end
    end
    % Determine the final minimum depth
    if (j<=length(d) && dmin_tur_num(i)<d(j))
        dmin_tur_num(i) = d(j);   
    end     
    if (j>length(d))
        dmin_tur_num(i) = NaN;
    end   
end

% filter out (d, u) where d_T>100 k (should not happen in a turbulent regime)
[row, col] = find(d_T>100);                         % record index of delta T where it is unreasonably large
for i=1:length(row)
    if (dmin_tur_num(row(i))<d(col(i)))
        dmin_tur_num(row(i))=d(col(i))+d(2)-d(1);
    end
end

% Disregard dmax and dmin_tur_num where dmax < dmin_tur_num
dmaxnum0 = dmax_num; 
dmax_num(dmax_num < dmin_tur_num) = NaN;
dmin_tur_num((dmaxnum0 < dmin_tur_num)|(isnan(dmaxnum0))) = NaN;

% Extract delta T and delta S at dmin
dT = zeros(1,length(u));
for i=1:length(u)
    if (ind_dmin(i))
        dT(i) = d_T(i, ind_dmin(i));
        s2_atmin = s2(i, ind_dmin(i));
    else
        dT(i) = nan;
        s2_atmin = nan;
    end
end

end