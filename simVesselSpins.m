function [magMap,vMap,pVessel,pSim,pMri] = simVesselSpins(pVessel, pSim, pMri, posFE, posPE)
if ~exist('posFE','var') || isempty(posFE); posFE = 0; end
if ~exist('posPE','var') || isempty(posPE); posPE = 0; end

% Define radial coordinates (relative to vessel center)
[gridFE, gridPE] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);  % dim1=FE, dim2=PE
rGrid = sqrt((gridFE - posFE).^2 + (gridPE - posPE).^2);

% Define compartments masks
% if ~isfield(pVessel,'mask'); pVessel.mask = struct('lumen',[],'wall',[],'surround',[]); end
% if isempty(pVessel.mask.lumen) || isempty(pVessel.mask.wall) || isempty(pVessel.mask.surround)
    pVessel.mask.lumen      = rGrid<=(pVessel.ID/2);                                    % vessel lumen
    pVessel.mask.wall       = rGrid> (pVessel.ID/2) & rGrid<=(pVessel.ID/2+pVessel.WT); % vessel wall
    pVessel.mask.surround   = rGrid> (pVessel.ID/2+pVessel.WT);                         % static surround
% end
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
% map signal magnitude
magMap = zeros(size(rGrid));
magMap(pVessel.mask.lumen)    = pVessel.S.lumen;
magMap(pVessel.mask.wall)     = pVessel.S.wall;
magMap(pVessel.mask.surround) = pVessel.S.surround;
nSpinPerVox = pSim.nSpinPerVox;
magMap = magMap ./ nSpinPerVox; % divide so summing spins in center voxel gives the measured signal


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

