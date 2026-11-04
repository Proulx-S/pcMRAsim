function results = test_inflowCount(verbose)
% test_inflowCount  Validate the per-isochromat inflow model of simVesselSpins (3D slab, RF pulse counting).
%
%   results = test_inflowCount        % prints a table, errors if a check fails
%   results = test_inflowCount(false) % quiet
%
% Needs util (getMxy_ss) on the path.
%
% Reference ("truth"): spins flow at constant velocity v through a rectangular slab, entering fully
% relaxed through its face; pulses every TR. Each spin's earlier pulses are counted from its entry
% time and its Mz is evolved explicitly, pulse by pulse. The column average of Mz just before a
% pulse is what a voxel spanning the slab measures.
%
% Checks (target protocol: 7 T 2D spoiled GRE, TR 8.4 ms, FA 40, slice 1.2 mm, voxel 0.4x0.4x1.2 mm)
%   1. Central isochromat column of a perpendicular vessel = pulse-by-pulse column average.
%   2. The result does not depend on the number of isochromat planes through the slab.
%   3. A vessel axis pointing to -z gives the mirrored column (same mean), and Mz <= M0.
%   4. A tilted vessel depends only on the through-plane velocity v*n_hat(3).
%   5. The surround (S.surround empty) uses the tissue T1.
%   6. getMz_n refuses non-integer pulse counts, and is unchanged for integer ones.
%   7. No NaN in the signal map.
% See docs/inflowPulseCount.md.

if ~exist('verbose','var') || isempty(verbose); verbose = true; end
assert(exist('getMxy_ss','file') == 2, 'test_inflowCount:util', 'util (getMxy_ss) must be on the path');
tolTruth = 2e-5;   % [M0] absolute; the reference's stratified-sampling error

[pMri, vCrit] = defaultMri();
k   = sind(pMri.FA)*exp(-pMri.TE/pMri.relax.blood.T2star);   % Mxy/Mz for blood
E1  = exp(-pMri.TR/pMri.relax.blood.T1);  Q1 = E1*cosd(pMri.FA);
Mss = (1-E1)/(1-Q1);

% --- 1. central column vs pulse-by-pulse reference (default grid: 9 planes) -----------------
frac = [0.03 0.05 0.1 0.2 0.3 0.45 0.5 0.7 0.9 0.99 1 1.5 2];
col  = nan(size(frac));  truth = col;  old = col;
for i = 1:numel(frac)
    [mMap, info] = runOne(pMri, 0.8, [0 0 1], frac(i)*vCrit/2, 2);   % parabolic1: centre v = 2*vMean
    assert(all(isfinite(mMap(:))), 'test_inflowCount:nan', 'NaN in mMap');              % check 7
    col(i)   = mean(mMap(info.ic,info.jc,:))/k;
    truth(i) = columnTruth(frac(i)*vCrit, pMri);
    old(i)   = oldColumn(frac(i)*vCrit, pMri, info.coorSLC);
end
assert(all(abs(col - truth) < tolTruth), 'test_inflowCount:truth', ...
    'central column differs from the pulse-by-pulse reference (max %.2g M0)', max(abs(col - truth)));

% --- 2. plane-count independence ----------------------------------------------------------------
isoMin = [1 2 4 7];  colN = nan(numel(isoMin), 3);  nPlanes = nan(size(isoMin));
fracN  = [0.05 0.3 0.9];
for a = 1:numel(isoMin)
    for b = 1:numel(fracN)
        [mMap, info] = runOne(pMri, 0.8, [0 0 1], fracN(b)*vCrit/2, isoMin(a));
        colN(a,b)  = mean(mMap(info.ic,info.jc,:))/k;
    end
    nPlanes(a) = numel(info.coorSLC);
end
spreadN = max(max(colN) - min(colN));
assert(spreadN < 1e-12, 'test_inflowCount:planes', 'result depends on the plane count (spread %.2g)', spreadN);

% --- 3. axis pointing to -z ---------------------------------------------------------------------
[mUp, info] = runOne(pMri, 0.8, [0 0  1], 0.3*vCrit/2, 2);
mDn         = runOne(pMri, 0.8, [0 0 -1], 0.3*vCrit/2, 2);
cUp = squeeze(mUp(info.ic,info.jc,:));  cDn = squeeze(mDn(info.ic,info.jc,:));
assert(max(abs(cUp - flipud(cDn))) < 1e-12, 'test_inflowCount:sign', '-z axis is not the mirrored column');
assert(max(mDn(:))/k <= 1 + 1e-12, 'test_inflowCount:signM0', '-z axis gives Mz > M0');

% --- 4. tilt: only v*n_hat(3) matters (wide vessel, so v ~ Vmax over the FOV) --------------------
th  = 30;  vC = 0.3*vCrit;
mT  = runOne(pMri, 200, [sind(th) 0 cosd(th)], vC/cosd(th)/2, 2);
mP  = runOne(pMri, 200, [0 0 1],                vC/2,          2);
cT  = mean(mT(info.ic,info.jc,:))/k;  cP = mean(mP(info.ic,info.jc,:))/k;
assert(abs(cT - cP) < 1e-4, 'test_inflowCount:tilt', 'tilted vessel differs from perpendicular at equal v*n_hat(3) (%.2g)', abs(cT-cP));

% --- 5. surround uses the tissue T1 -------------------------------------------------------------
[mS, ~, tissueMz] = runOne(pMri, 0.8, [0 0 1], 0.3*vCrit/2, 2, []);
E1g = exp(-pMri.TR/pMri.relax.GM.T1);
MssGM = (1-E1g)/(1-E1g*cosd(pMri.FA));
assert(abs(tissueMz - MssGM) < 1e-12, 'test_inflowCount:surround', 'surround Mz %.5f, expected GM steady state %.5f', tissueMz, MssGM);
assert(all(isfinite(mS(:))), 'test_inflowCount:nanS', 'NaN in mMap');

% --- 6. getMz_n: integer counts only, unchanged for integers --------------------------------
n = [0 1 2 5 40 inf];
Mz = getMz_n(pMri, pMri.relax.blood, n);
assert(max(abs(Mz - (Mss + (1-Mss)*Q1.^n))) < 1e-14, 'test_inflowCount:integerN', 'getMz_n changed for integer n');
threw = false;
try
    getMz_n(pMri, pMri.relax.blood, 2.5);
catch ME
    threw = strcmp(ME.identifier,'getMz_n:nonIntegerN');
end
assert(threw, 'test_inflowCount:nonIntegerN', 'getMz_n accepted a non-integer pulse count');

% --- report --------------------------------------------------------------------------------------
results = struct('frac',frac,'truth',truth,'new',col,'old',old, 'isoMin',isoMin,'nPlanes',nPlanes, ...
    'colN',colN,'spreadN',spreadN,'tilt',[cT cP],'surroundMz',tissueMz);
if verbose
    fprintf('\nCentral column of a perpendicular vessel, Mz/M0 (7 T, TR 8.4 ms, FA 40, 1.2 mm slab, 9 planes)\n');
    fprintf('  Q1 = %.4f, Mz(v=0) = %.4f, vCrit = %.2f cm/s\n', Q1, Mss, vCrit);
    fprintf('  v/vCrit   truth    new      |err|        previous (continuous n)   error of the excess\n');
    for i = 1:numel(frac)
        fprintf('  %6.2f   %7.4f  %7.4f  %8.1e      %7.4f                  %+6.1f %%\n', frac(i), truth(i), col(i), ...
            abs(col(i)-truth(i)), old(i), 100*(old(i)-truth(i))/max(truth(i)-Mss,eps));
    end
    fprintf('\nPlane-count independence: planes through the slab = %s; max spread of the column mean = %.1e\n', ...
        mat2str(nPlanes), spreadN);
    fprintf('Axis -z: mirrored column (max diff %.1e), max Mz/M0 = %.4f\n', max(abs(cUp - flipud(cDn))), max(mDn(:))/k);
    fprintf('Tilt %d deg at equal v*n_hat(3): %.6f vs perpendicular %.6f\n', th, cT, cP);
    fprintf('Surround Mz = %.5f (GM steady state %.5f; blood would give %.5f)\n', tissueMz, MssGM, Mss);
    fprintf('\ntest_inflowCount: all checks passed\n');
end
end


function [pMri, vCrit] = defaultMri()
pMri.TR = 8.4e-3; pMri.TE = 3.68e-3; pMri.FA = 40; pMri.sliceThickness = 1.2;
pMri.relax.blood.T1 = 2.58;  pMri.relax.blood.T2star = 10e-3;
pMri.relax.GM.T1    = 1.939; pMri.relax.GM.T2star    = 32.9e-3;
vCrit = 0.1*pMri.sliceThickness/pMri.TR;   % [cm/s]
end


function [mMap, info, tissueMz] = runOne(pMri, ID, n_hat, vMean, isoMin, surround)
% 3x3 voxels of 0.4x0.4x1.2 mm, vessel through the centre
if ~exist('surround','var'); surround = 0; end
fov = [1.2 1.2 pMri.sliceThickness]; mat = [3 3 1];
[pSim.voxGrid, pSim.spinGrid, pSim.isochromatPerVox] = setGrid(fov, mat, isoMin, 'pseudoVoxel');
pVessel = struct('ID',ID,'WT',0,'profile','parabolic1','pos',[0 0 0],'n_hat',n_hat./norm(n_hat),'vMean',vMean);
pVessel.S.surround = surround; pVessel.S.wall = 0;
mMap = simVesselSpins(pVessel, pSim, pMri);
info.ic = find(abs(pSim.spinGrid.coorFE) < 1e-9);
info.jc = find(abs(pSim.spinGrid.coorPE) < 1e-9);
info.coorSLC = pSim.spinGrid.coorSLC;
tissueMz = mMap(1,1,1)/(sind(pMri.FA)*exp(-pMri.TE/pMri.relax.GM.T2star));   % corner: tissue
end


function Mbar = columnTruth(v, pMri)
% Pulse-by-pulse column average of Mz just before a pulse (M0 = 1), plug flow at v [cm/s].
nSpin = 200003;
E1 = exp(-pMri.TR/pMri.relax.blood.T1);
tTransit = pMri.sliceThickness/(10*v);
nPulse   = ceil(tTransit/pMri.TR) + 3;
tObs     = nPulse*pMri.TR;
tEntry   = tObs - ((1:nSpin) - 0.5)/nSpin*tTransit;
nPrior   = nPulse - floor(tEntry/pMri.TR) - 1;
Mz = ones(1, nSpin);
for k = 1:max(nPrior)
    hit = nPrior >= k;
    Mz(hit) = Mz(hit)*cosd(pMri.FA)*E1 + (1 - E1);
end
Mbar = mean(Mz);
end


function Mbar = oldColumn(v, pMri, coorSLC)
% The previous implementation, for comparison: depth from the first plane, continuous n in Q1^n.
E1 = exp(-pMri.TR/pMri.relax.blood.T1); Q1 = E1*cosd(pMri.FA); Mss = (1-E1)/(1-Q1);
n  = (coorSLC - coorSLC(1))/10/v/pMri.TR;
Mbar = mean(Mss + (1-Mss)*Q1.^n);
end
