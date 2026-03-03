function [magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, x0, y0)
if ~exist('x0','var') || isempty(x0); x0 = pVessel.x0; end
if ~exist('y0','var') || isempty(y0); y0 = pVessel.y0; end


% Define vessel radial coordinates
rGrid = sqrt((xGrid-x0).^2 + (yGrid-y0).^2);
nSpin = numel(rGrid);


% Define compartments masks
mask.lumen      = rGrid<=(pVessel.ID/2); % vessel lumen
mask.wall       = rGrid> (pVessel.ID/2) & rGrid<=(pVessel.ID/2+pVessel.WT); % vessel wall
mask.surround   = rGrid> (pVessel.ID/2+pVessel.WT); % static surround
mask.lumenPlug  = rGrid<(pVessel.PD/2) & mask.lumen; % plug flow center
mask.lumenLami  = rGrid>=(pVessel.PD/2) & mask.lumen; % laminar flow region


% Define ROI compartment fractions
cmptList = fieldnames(mask);
for iCmpt = 1:length(cmptList)
    pVessel.volFraction.(cmptList{iCmpt}) = mean(mask.(cmptList{iCmpt}),[1 2]);
end


% Define spin velocity map
vMap = getVelMap(rGrid, pVessel.ID, pVessel.profile, pVessel.PD); % [cm/s]
if ~isempty(pVessel.vMax) && isempty(pVessel.vMean)
    vMap = scale2maxVel(vMap, pVessel.vMax); % to the desired maximum velocity
elseif ~isempty(pVessel.vMean) && isempty(pVessel.vMax)
    vMap = scale2meanVel(vMap, pVessel.vMean, mask.lumen); % to the desired mean velocity
else
    error('Either pVessel.vMax or pVessel.vMean must be specified');
end


% Define spin signal magnitude
cmptList = fieldnames(pVessel.cmptMag);
for iCmpt = 1:length(cmptList)
    pVessel.spinMag.(cmptList{iCmpt}) = pVessel.cmptMag.(cmptList{iCmpt})./nSpin;
end

% Define signal magnitude map
magMap = zeros(size(rGrid));
if length(pVessel.spinMag.lumenLami) == 1
    magMap(mask.lumenLami) = pVessel.spinMag.lumenLami; % laminar flow region of the vessel lumen (single spin magnitudes)
else
    %                   vq = interp1(x                             ,v                             ,xq  )
    magMap(mask.lumenLami) = interp1(pVessel.cmptMag.lumenLami(1,:),pVessel.cmptMag.lumenLami(2,:),vMap(mask.lumenLami),'linear','extrap');
end
magMap(mask.lumenPlug) = pVessel.spinMag.lumenPlug; % plug flow center of the vessel lumen (single spin magnitudes)
magMap(mask.wall)      = pVessel.spinMag.wall; % wall of the vessel (single spin magnitudes)
magMap(mask.surround)  = pVessel.spinMag.surround; % static surround of the vessel (single spin magnitudes)

