function [magMap,vMap,pVessel,pSim,pMri] = simVesselSpins(pVessel, pSim, pMri, posFE, posPE)
if ~exist('posFE','var') || isempty(posFE); posFE = 0; end
if ~exist('posPE','var') || isempty(posPE); posPE = 0; end

% Define radial coordinates (relative to vessel center)
rGrid = sqrt((pSim.gridFE-posFE).^2 + (pSim.gridPE-posPE).^2);

% Define compartments masks
% if ~isfield(pVessel,'mask'); pVessel.mask = struct('lumen',[],'wall',[],'surround',[]); end
% if isempty(pVessel.mask.lumen) || isempty(pVessel.mask.wall) || isempty(pVessel.mask.surround)
    pVessel.mask.lumen      = rGrid<=(pVessel.ID/2);                                    % vessel lumen
    pVessel.mask.wall       = rGrid> (pVessel.ID/2) & rGrid<=(pVessel.ID/2+pVessel.WT); % vessel wall
    pVessel.mask.surround   = rGrid> (pVessel.ID/2+pVessel.WT);                         % static surround
% end

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
% vessel surround (static)
Mxy = getMxy_ss(getMz_ss(pMri,pMri.relax.GM),pMri,pMri.relax.GM);
pVessel.S.surround = Mxy;
% map signal magnitude
magMap = zeros(size(rGrid));
magMap(pVessel.mask.lumen)    = pVessel.S.lumen;
magMap(pVessel.mask.surround) = pVessel.S.surround;
magMap = magMap./pSim.nSpin; % divide by the number of spins a voxel, so summing the spins gives the measured signal in a voxel


% Precompute montecarlo tessalation
if pSim.monteCarloN > 0 && (~isfield(pSim,'monteCarloShiftFE') || ~isfield(pSim,'monteCarloShiftPE') || isempty(pSim.monteCarloShiftFE) || isempty(pSim.monteCarloShiftPE))
    nSpinFE = max(sum(pSim.gridVoxIdx==0,2));
    shiftFE = (1:nSpinFE)-nSpinFE/2-0.5;
    nSpinPE = max(sum(pSim.gridVoxIdx==0,1));
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

