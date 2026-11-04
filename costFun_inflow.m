function [residual,sim] = costFun_inflow(theta, pVessel, pSim, pMri, data)
% costFun_inflow  Residuals for the cylinder3D inflow forward model.
%
% theta = [Vmax, R, nx, ny, cx_FE, cx_PE, cx_SLC, A, WT, S_tissue, sigma_n]
%
% data fields:
%   pSim_base    - pSim with voxGrid.fovFE/PE, matFE/PE, nSpin, monteCarloN set
%   pMri_base    - pMri with TR, TE, FA, sliceThickness, relax.blood
%   m_meas       - [nFE x nPE] measured magnitude (venc=inf average)
%   v_meas       - [nFE x nPE] measured velocity (best venc)
%   mask_vel     - [nFE x nPE] logical mask for velocity residuals (blood pixels)
%   m_noflow     - [nFE x nPE] noFlow magnitude, or [] to skip noFlow residuals
%   noise_grid   - [nFE x nPE] complex noise realization for flow data
%   noise_noflow - [nFE x nPE] complex noise realization for noFlow data
%   FEgrid       - [nFE x nPE] FE coordinates of voxel centers [mm]
%   PEgrid       - [nFE x nPE] PE coordinates of voxel centers [mm]
%   sv           - velocity residual scale factor [cm/s]
%   sm           - magnitude residual scale factor [a.u.]

Vmax     = theta(1); %
R        = theta(2); %
nx       = theta(3); %
ny       = theta(4); %
cx_FE    = theta(5); %
cx_PE    = theta(6); %
cx_SLC   = theta(7); %
A        = theta(8); %
WT       = theta(9); %
S_tissue = theta(10); %
sigma_n  = theta(11); %
nz = sqrt(max(0, 1 - nx^2 - ny^2)); % ||||Note for Claude: we should probably reparameterize this in degrees tilt to avoid imposible vectors||||

% --- inject in sim structure ---
pVessel.ID         = 2 * R;
pVessel.WT         = WT;
% Erratum pcMRAsim#4: the simulation reads pVessel.pos and pVessel.vMean. This used to set
% posFE/posPE/posSLC and Vmax, which nothing reads, so the fit was flat in position and velocity.
pVessel.pos        = [cx_FE cx_PE cx_SLC];
pVessel.n_hat      = [nx; ny; nz];
switch pVessel.profile
    case 'parabolic1'
        pVessel.vMean = Vmax/2;          % simVesselSpins: Vmax = 2*vMean
    otherwise
        error('costFun_inflow:profile', 'Vmax -> vMean is only defined for profile ''parabolic1'' (got ''%s'').', pVessel.profile);
end
pVessel.S.surround = S_tissue;
% A (theta(8)) and sigma_n (theta(11)) are NOT used by the forward model: no live code reads an
% amplitude scale (pSim.mriScale is only a default) or a noise scale. Their fitted values are
% meaningless until the model implements them (pcMRAsim#4).
persistent warnedInert
if isempty(warnedInert)
    warning('costFun_inflow:inertParams', ['theta(8) (A) and theta(11) (sigma_n) do not affect the forward model; ' ...
        'their fitted values are meaningless (pcMRAsim#4).']);
    warnedInert = true;
end
pVessel.A          = A;                  % kept for the caller's bookkeeping only
pSim.spinGrid.gridNoiseScale = sigma_n;  % idem


sim = runSim(pVessel, pSim, pMri, false, true, false);

residual = data(:)-sim.I(:);
residual = [real(residual); imag(residual)];






% % --- Aggregate per-spin magMap → per-voxel Mxy ---
% % magMap is [nTotalFE, nTotalPE] = [sFE*nVxFE, sPE*nVxPE].
% % Spins within each voxel are contiguous in each direction (setGrid uses uniform spacing).
% nVxFE = numel(res_sim.pSim.voxGrid.coorFE);
% nVxPE = numel(res_sim.pSim.voxGrid.coorPE);
% sFE   = round(numel(res_sim.pSim.spinGrid.coorFE) / nVxFE);
% sPE   = round(numel(res_sim.pSim.spinGrid.coorPE) / nVxPE);
% Mxy   = squeeze(sum(sum(reshape(double(res_sim.magMap), [sFE, nVxFE, sPE, nVxPE]), 1), 3));

% % --- Predicted magnitude with pre-realized complex noise ---
% m_pred = sqrt((Mxy + sigma_n .* real(data.noise_grid)).^2 + ...
%               (sigma_n .* imag(data.noise_grid)).^2);

% % --- Magnitude residuals (all voxels) ---
% res_mag = (data.m_meas(:) - m_pred(:)) / data.sm;

% % --- Velocity residuals at voxel centers (Z = 0, slab midplane) ---
% fe_rel  = data.FEgrid - cx_FE;
% pe_rel  = data.PEgrid - cx_PE;
% z_rel   = 0 - cx_SLC;
% dot_n   = fe_rel.*nx + pe_rel.*ny + z_rel.*nz;
% r_perp2 = fe_rel.^2 + pe_rel.^2 + z_rel.^2 - dot_n.^2;
% v_pred  = Vmax .* max(0, 1 - r_perp2 ./ max(R^2, eps));
% lumen_px = r_perp2 <= R^2;
% mask_v   = data.mask_vel & lumen_px;
% res_vel  = [];
% if any(mask_v(:))
%     res_vel = (data.v_meas(mask_v) - v_pred(mask_v)) / data.sv;
% end

% % --- noFlow magnitude residuals ---
% res_nf = [];
% if ~isempty(data.m_noflow)
%     pVessel_nf      = pVessel;
%     pVessel_nf.Vmax = 0;
%     res_nf_sim = runSim(pVessel_nf, data.pSim_base, data.pMri_base, false, false, true);
%     Mxy_nf = squeeze(sum(sum(reshape(double(res_nf_sim.magMap), [sFE, nVxFE, sPE, nVxPE]), 1), 3));
%     m_pred_nf = sqrt((Mxy_nf + sigma_n .* real(data.noise_noflow)).^2 + ...
%                      (sigma_n .* imag(data.noise_noflow)).^2);
%     res_nf = (data.m_noflow(:) - m_pred_nf(:)) / data.sm;
% end

% res = double([res_mag; res_vel(:); res_nf(:)]);
% end
