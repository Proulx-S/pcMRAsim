function [magMap,vMap,mask,pVessel] = simVesselSpins(pVessel, pSim, pMri,posFE,posPE)
if ~exist('posFE','var') || isempty(posFE); posFE = pVessel.posFE; end
if ~exist('posPE','var') || isempty(posPE); posPE = pVessel.posPE; end

% Define radial coordinates (relative to vessel center)
rGrid = sqrt((pSim.gridFE-posFE).^2 + (pSim.gridPE-posPE).^2);

% Define compartments masks
mask.lumen      = rGrid<=(pVessel.ID/2); % vessel lumen
mask.wall       = rGrid> (pVessel.ID/2) & rGrid<=(pVessel.ID/2+pVessel.WT); % vessel wall
mask.surround   = rGrid> (pVessel.ID/2+pVessel.WT); % static surround

% Define spin velocity map
vMap = getVelMap(rGrid, pVessel.ID, pVessel.profile, pVessel.PD); % [cm/s]
if ~isempty(pVessel.vMax) && isempty(pVessel.vMean)
    vMap = scale2maxVel(vMap, pVessel.vMax); % to the desired maximum velocity
elseif ~isempty(pVessel.vMean) && isempty(pVessel.vMax)
    vMap = scale2meanVel(vMap, pVessel.vMean, mask.lumen); % to the desired mean velocity
else
    error('Either pVessel.vMax or pVessel.vMean must be specified');
end


% MR signal magnitude
if isempty(pVessel.S.lumen)
    inflow.sliceThickness = pMri.sliceThickness;
    inflow.FA             = pMri.FA;
    inflow.T1             = pMri.relax.blood.T1;
    inflow.TR             = pMri.TR;
    inflow.T2star = pMri.relax.blood.T2star;
    inflow.TE     = pMri.TE;
    inflow.vMean  = pVessel.vMean;
    if isempty(inflow.vMean) && ~isempty(pVessel.vMax) && strcmp(pVessel.profile,'parabolic1')
        inflow.vMean   = pVessel.vMax/2;
    end
    [inflow.fAtvMean, inflow.Mz_v0, inflow.fMax, inflow.vCrit] = inflowEnhancementBianciardi(inflow.vMean,inflow);
    

    dbstack;
    keyboard;
    %!!!!!!!!! Implement transverse relaxation into inflow enhancement




    inflow.vCrit  = pVessel.vCrit;
    inflow.fMax   = fMax;
    pVessel.S.lumen
end


% !!!!!!!!!!!!
S.lumen = pVessel.S.lumen;
S.wall  = pVessel.S.wall;
S.surround = pVessel.S.surround;

% Define signal magnitude map
magMap = zeros(size(rGrid));
f      = nan(size(vMap));
vel    = vMap(mask.lumen);
pVessel
[f(mask.lumen), Mz_v0, fMax, vCrit] = inflowEnhancementBianciardi(vel,pVessel.inflow);
Mxy_v0 = Mz_v0.*sin(pSim.FA*pi/180).*exp(-pSim.TE./pVessel.T2star.blood);
magMap(mask.lumen)      = f(mask.lumen).*Mxy_v0;
pVessel.inflow.Mz_v0    = Mz_v0;
pVessel.inflow.fMax     = fMax;
pVessel.inflow.vCrit    = vCrit;

if nnz(mask.lumenPlug)
    dbstack; error('double-check that');
    magMap(mask.lumenPlug) = pVessel.S.lumenPlug; % plug flow center of the vessel lumen (single spin magnitudes)
end
if nnz(mask.wall)
    dbstack; error('double-check that');
    magMap(mask.wall) = pVessel.S.wall; % wall of the vessel (single spin magnitudes)
end
magMap(mask.surround) = pVessel.S.surround; % static surround of the vessel (single spin magnitudes)
magMap = magMap./nSpin; % MR signal between 0 and 1 for the full ROI (not voxel). Dividing by the number of spins in the ROI for the signal of individual spins.

