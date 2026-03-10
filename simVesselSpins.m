function [magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, pSim, x0, y0, anaFlag, vFix)
if ~exist('x0','var')      || isempty(x0); x0 = pVessel.x0; end
if ~exist('y0','var')      || isempty(y0); y0 = pVessel.y0; end
if ~exist('anaFlag','var') || isempty(anaFlag); anaFlag = 'inflowOnSpinVelocity'; end
    % anaFlag = 'inflowFixedAtMax';
    % anaFlag = 'inflowFixedAtVelocity';
    % anaFlag = 'inflowOnMeanVelocity';
    % anaFlag = 'inflowOnSpinVelocity';

% NOTE: The 0 to 1 MR signal is defined on the full ROI
%       Single-spin signal is the full ROI signal divided by the number of spins in the ROI
%       This part of the code is not aware of voxels.
%       This signal definition needs to be taken into account when later summing signal in a voxel (ROI subset). 


% Define vessel radial coordinates
rGrid = sqrt((xGrid-x0).^2 + (yGrid-y0).^2);
nSpin = numel(rGrid);


% Define compartments masks
mask.lumen      = rGrid<=(pVessel.ID/2); % vessel lumen
mask.wall       = rGrid> (pVessel.ID/2) & rGrid<=(pVessel.ID/2+pVessel.WT); % vessel wall
mask.surround   = rGrid> (pVessel.ID/2+pVessel.WT); % static surround
mask.lumenPlug  = rGrid< (pVessel.PD/2) & mask.lumen; % plug flow center
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


% % Define spin signal magnitude
% cmptList = fieldnames(pVessel.S);
% for iCmpt = 1:length(cmptList)
%     pVessel.spinMag.(cmptList{iCmpt}) = pVessel.S.(cmptList{iCmpt});
% end

% Define signal magnitude map
magMap = zeros(size(rGrid));
f      = nan(size(vMap));
switch anaFlag
    case 'inflowFixedAtMax'
        vel = inf;
    case 'inflowFixedAtVelocity'
        vel = vFix;
    case 'inflowOnMeanVelocity'
        vel = pVessel.vMean;
    case 'inflowOnSpinVelocity'
        vel = vMap(mask.lumen);
    otherwise
        dbstack; error('invalid anaFlag');
end
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

