function [magMap,vMap,pVessel,pSim,pMri] = simVesselSpins(pVessel, pSim, pMri)

% Define radial coordinates (relative to vessel center)
[gridFE, gridPE] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);  % dim1=FE, dim2=PE
rGrid = sqrt((gridFE - pVessel.posFE).^2 + (gridPE - pVessel.posPE).^2);

% Define compartment masks
% Elliptical vessel if pVessel.AR and pVessel.alpha are provided; circular otherwise.
R_lumen = pVessel.ID / 2;  % semi-major axis [mm]
if isfield(pVessel,'AR') && ~isempty(pVessel.AR) && pVessel.AR ~= 1
    alpha_v = pVessel.alpha;
    AR_v    = pVessel.AR;
    dFE = gridFE - pVessel.posFE;
    dPE = gridPE - pVessel.posPE;
    u   =  dPE.*cos(alpha_v) + dFE.*sin(alpha_v);   % along major axis
    w   = -dPE.*sin(alpha_v) + dFE.*cos(alpha_v);   % along minor axis
    ellipseDist = sqrt(u.^2 + (AR_v .* w).^2);      % = R_lumen at ellipse boundary
    pVessel.mask.lumen    = ellipseDist <= R_lumen;
    pVessel.mask.wall     = ellipseDist >  R_lumen & ellipseDist <= R_lumen + pVessel.WT;
    pVessel.mask.surround = ~pVessel.mask.lumen & ~pVessel.mask.wall;
else
    pVessel.mask.lumen    = rGrid <= R_lumen;
    pVessel.mask.wall     = rGrid >  R_lumen & rGrid <= R_lumen + pVessel.WT;
    pVessel.mask.surround = rGrid >  R_lumen + pVessel.WT;
end
assert(all(pVessel.mask.lumen(:) + pVessel.mask.wall(:) + pVessel.mask.surround(:) == 1), ...
    'pVessel.mask: lumen, wall, surround must be non-overlapping and cover all spins');

% Define spin velocity map
if ischar(pVessel.profile)
    vMap = getVelMap(rGrid, pVessel.ID, pVessel.profile, pVessel.PD); % [cm/s]
    if ~isempty(pVessel.vMax) && isempty(pVessel.vMean)
        vMap = scale2maxVel(vMap, pVessel.vMax); % to the desired maximum velocity
    elseif ~isempty(pVessel.vMean) && isempty(pVessel.vMax)
        vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
    elseif strcmp(pVessel.profile,'parabolic1') && ~isempty(pVessel.vMax) && ~isempty(pVessel.vMean) && pVessel.vMax/2==pVessel.vMean
        vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
    else
        error('Either pVessel.vMax or pVessel.vMean must be specified');
    end
elseif isnumeric(pVessel.profile)
    vMap = zeros(size(rGrid));
    vMap(:) = pVessel.profile;
else
    dbstack; error('Invalid vessel profile');
end


% MR signal magnitude
% vessel lumen signal (flowing)
if isempty(pVessel.S.lumen)
    switch pVessel.profile
        case 'plug'
            [Mz_vMean,pMri] = getMz_ss(          pMri,pMri.relax.blood,pVessel.vMean);
            [Mxy_vMax,pMri] = getMxy_ss(Mz_vMean,pMri,pMri.relax.blood              );
            pVessel.S.lumen = Mxy_vMax;
        case {'parabolic','parabolic1'}
            [Mz ,pMri] = getMz_ss(    pMri,pMri.relax.blood,vMap(pVessel.mask.lumen));
            [Mxy,pMri] = getMxy_ss(Mz,pMri,pMri.relax.blood                         );
            pVessel.S.lumen = Mxy;
        otherwise
            dbstack; error('Invalid vessel profile');
    end
end
% vessel surround (static)
if isempty(pVessel.S.surround)
    pVessel.S.surround = getMxy_ss(getMz_ss(pMri,pMri.relax.GM),pMri,pMri.relax.GM);
end
% map signal magnitude (divide by nSpinPerVox so summing spins within voxels gives S, the hypothetical measured signal if the voxel was single-compartment)
magMap = zeros(size(rGrid));
magMap(pVessel.mask.lumen)    = pVessel.S.lumen    ./pSim.nSpinPerVox;
magMap(pVessel.mask.wall)     = pVessel.S.wall     ./pSim.nSpinPerVox;
magMap(pVessel.mask.surround) = pVessel.S.surround ./pSim.nSpinPerVox;


% Precompute montecarlo tessalation
if pSim.monteCarloN > 0 && (~isfield(pSim,'monteCarloShiftFE') || ~isfield(pSim,'monteCarloShiftPE') || isempty(pSim.monteCarloShiftFE) || isempty(pSim.monteCarloShiftPE))
    nSpinFE = pSim.spinGrid.matFE / pSim.voxGrid.matFE;  % spins per voxel in FE
    shiftFE = (1:nSpinFE)-nSpinFE/2-0.5;
    nSpinPE = pSim.spinGrid.matPE / pSim.voxGrid.matPE;  % spins per voxel in PE
    shiftPE = (1:nSpinPE)-nSpinPE/2-0.5;
    % find all possible combination of FE and PE shifts
    [idx1, idx2] = ndgrid(1:length(shiftFE), 1:length(shiftPE));
    idx = [idx1(:), idx2(:)];
    % remove the no-shift combination since it is always done before
    idx(all(idx==[find(shiftFE==0) find(shiftPE==0)],2),:) = [];
    % shuffle
    if pSim.monteCarloN~=inf
        idx = idx(randperm(size(idx,1),pSim.monteCarloN),:);
    end
    pSim.monteCarloShiftFE = shiftFE(idx(:,1));
    pSim.monteCarloShiftPE = shiftPE(idx(:,2));
end

