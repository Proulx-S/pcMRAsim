function test_costFun_inflow
% test_costFun_inflow  Every parameter costFun_inflow passes to the forward model must move the
% residual (pcMRAsim#4: position and velocity used to be written to fields nothing read).
%
%   addpath(pcMRAsimDir); addpath(fullfile(pcMRAsimDir,'tests')); addpath(utilDir); test_costFun_inflow

p = runSim;                                   % defaults: 3 x 3 voxels, 7 T
pVessel = p.pVessel; pSim = p.pSim; pMri = p.pMri;
pSim.gridMode = 'allVox';                     % as multiVencISMRM2026 doIt_inflow section 13
pSim.voxGrid.fovSLC = pMri.sliceThickness;    % idem (older runSim defaults differ: 2 vs 1 mm)
pVessel.geometry = 'infiniteCylinder';
pVessel.ID = 0.6; pVessel.WT = 0.05; pVessel.vMean = 2;
pVessel.n_hat = [0; 0; 1];
pVessel.S.surround = 0.02;
if isempty(pVessel.S.wall); pVessel.S.wall = 0; end

%        [Vmax, R  , nx, ny, cx_FE, cx_PE, cx_SLC, A, WT  , S_tissue, sigma_n]
theta1 = [ 4  , 0.3, 0 , 0 , 0.10 , -0.05, 0     , 1, 0.05, 0.02    , 0      ];
[~, sim1] = costFun_inflow(theta1, pVessel, pSim, pMri, zeros(1));
data = sim1.I;
r0 = costFun_inflow(theta1, pVessel, pSim, pMri, data);
assert(max(abs(r0)) < 1e-12, 'residual at the generating theta is not zero');

names = {'Vmax','R','nx','ny','cx_FE','cx_PE','cx_SLC','A','WT','S_tissue','sigma_n'};
step  = [ 1    , 0.05, 0.05, 0.05, 0.1  , 0.07 , 0     , 0.5, 0.05, 0.01   , 0.1    ];
acts  = [ true , true, true, true, true , true , false , false, true, true , false  ];  % cx_SLC: slab is z-invariant for nz=1
% (cx_PE steps by 0.07, not 0.1: -0.05 -> +0.05 would be the mirror image on the 3x3 grid, an identical signal)
fprintf('%-9s %-12s %s\n', 'param', '|residual|', 'expected to act');
for k = 1:numel(names)
    if step(k) == 0; continue; end
    th = theta1; th(k) = th(k) + step(k);
    w = warning('off', 'costFun_inflow:inertParams');
    r = costFun_inflow(th, pVessel, pSim, pMri, data);
    warning(w);
    nr = norm(r);
    fprintf('%-9s %-12.3g %d\n', names{k}, nr, acts(k));
    if acts(k)
        assert(nr > 1e-9, '%s does not move the residual', names{k});
    else
        assert(nr < 1e-12, '%s unexpectedly moves the residual (update this test and pcMRAsim#4)', names{k});
    end
end
disp('test_costFun_inflow: all checks passed')
end
