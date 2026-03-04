function [magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, x0, y0, anaFlag)
if ~exist('x0','var')      || isempty(x0); x0 = pVessel.x0; end
if ~exist('y0','var')      || isempty(y0); y0 = pVessel.y0; end
if ~exist('anaFlag','var') || isempty(anaFlag); anaFlag = 'inflowOnSpinVelocity'; end
    % anaFlag = 'noInflow';
    % anaFlag = 'inflowOnMeanVelocity';
    % anaFlag = 'inflowOnSpinVelocity';


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
% cmptList = fieldnames(pVessel.cmptMag);
% for iCmpt = 1:length(cmptList)
%     pVessel.spinMag.(cmptList{iCmpt}) = pVessel.cmptMag.(cmptList{iCmpt});
% end

% Define signal magnitude map
magMap = zeros(size(rGrid));
switch anaFlag
    case 'noInflow'
        magMap(mask.lumenLami) = pVessel.cmptMag.lumenLami(end);
    case {'inflowOnMeanVelocity','inflowOnSpinVelocity'}
        f = nan(size(vMap));
        switch anaFlag
            case 'inflowOnMeanVelocity'
                [f(mask.lumenLami), Mz_v0, fMax, vCrit] = inflowEnhancementBianciardi(pVessel.vMean,pVessel.inflow);
            case 'inflowOnSpinVelocity'
                [f(mask.lumenLami), Mz_v0, fMax, vCrit] = inflowEnhancementBianciardi(vMap(mask.lumenLami),pVessel.inflow);
            otherwise
                dbstack; error('invalid anaFlag');
        end
        magMap(mask.lumenLami)  = f(mask.lumenLami).*Mz_v0;
        pVessel.inflow.Mz_v0    = Mz_v0;
        pVessel.inflow.fMax     = fMax;
        pVessel.inflow.vCrit    = vCrit;
    otherwise
        dbstack; error('invalid anaFlag');
end

if nnz(mask.lumenPlug)
    dbstack; error('double-check that');
    magMap(mask.lumenPlug) = pVessel.cmptMag.lumenPlug; % plug flow center of the vessel lumen (single spin magnitudes)
end
magMap(mask.wall)      = pVessel.cmptMag.wall; % wall of the vessel (single spin magnitudes)
magMap(mask.surround)  = pVessel.cmptMag.surround; % static surround of the vessel (single spin magnitudes)
magMap = magMap./nSpin;

