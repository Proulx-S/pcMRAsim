function [magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, pSim, x0, y0, anaFlag, vFix)
% simVesselSpins  Spin magnitude and velocity maps for a vessel cross-section, with inflow enhancement.
%
% INFLOW MODEL VERSION: continuous-u (Bianciardi et al. 2016 eq. 2.4 at non-integer u = vCrit/v).
%   The lumen inflow factor comes from inflowEnhancementBianciardi, which evaluates eq. 2.4,
%   (1 - q^u)/(u*(1 - q)), at u = vCrit/v for each spin's velocity (or the velocity anaFlag selects).
%   Eq. 2.4 is derived for integer u only; between integers this smooth curve is up to a few % high
%   compared with the exact slice-column average (piecewise linear in v). Not the recommended model:
%   the exact version is util getMz_ss on a util errata branch (util#15) and, per isochromat,
%   simVesselSpins/getMz_n on the pcMRAsim errata branch (pcMRAsim#5). Kept unchanged for
%   reproducibility of results computed with it.
%
% REFERENCES
%   Bianciardi M, Toschi N, Polimeni JR, et al. (2016) The pulsatility volume index: an indicator of
%     cerebrovascular compliance based on fast magnetic resonance imaging of cardiac and respiratory
%     pulsatility. Phil Trans R Soc A 374:20150184. doi:10.1098/rsta.2015.0184
%     (three-regime formulation; eq. 2.4, derived for integer vCrit/v)
%   Gullberg GT, Simons MA, Wehrli FW (1988) A mathematical model for signal from spins flowing during
%     the application of spin echo pulse sequences. Magn Reson Imaging 6:437-461.
%     doi:10.1016/0730-725X(88)90481-X (fractional-volume method)
%   Gao J-H, Holland SK, Gore JC (1988) Nuclear magnetic resonance signal from flowing nuclei in rapid
%     imaging using gradient echoes. Med Phys 15:809-814. doi:10.1118/1.596197
%     (that method for spoiled FLASH: exact for non-integer vCrit/v, sec. III.A, eq. 26)
%   Comparison of the four inflow versions (2026-11-06; private repo, the label above stands on its own):
%   https://github.com/Proulx-S/microscopyValidation/blob/9057997d65a96cad63419fed24daa1ba885ed6c6/docs/inflowBetweenIntegerPulseCounts.md
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

