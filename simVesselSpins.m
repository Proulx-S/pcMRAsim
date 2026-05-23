function [mMap,vMap,pVessel,pSim,pMri] = simVesselSpins(pVessel, pSim, pMri)


% --- isochromat coordinate grids ---
[gridFE, gridPE, gridSLC] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE, pSim.spinGrid.coorSLC);

% --- infinite cylinder radius grid ---
gridFE_vessel  = gridFE  - pVessel.pos(1);
gridPE_vessel  = gridPE  - pVessel.pos(2);
gridSLC_vessel = gridSLC - pVessel.pos(3);
dot_n  = gridFE_vessel.*pVessel.n_hat(1) + gridPE_vessel.*pVessel.n_hat(2) + gridSLC_vessel.*pVessel.n_hat(3);
gridR = sqrt(max(0, gridFE_vessel.^2 + gridPE_vessel.^2 + gridSLC_vessel.^2 - dot_n.^2));

% --- masks and vox belonging ---
lumenMask  =   gridR<=(pVessel.ID/2             )               ;
wallMask   = ( gridR< (pVessel.ID/2 + pVessel.WT) ) & ~lumenMask;
tissueMask =   gridR>=(pVessel.ID/2 + pVessel.WT)               ;
pSim.spinGrid.voxIdx = getVoxIdx(pSim.voxGrid, pSim.spinGrid, gridFE, gridPE, gridSLC, gridR);


% --- parabolic velocity along cylinder ---
switch pVessel.profile
    case 'parabolic1'
        Vmax  = pVessel.vMean*2;
    otherwise
        dbstack; error('double-check that')
end
vMap = zeros(size(gridR));
vMap(lumenMask) = Vmax .* max(0, (1 - (gridR(lumenMask)./(pVessel.ID/2)).^2) );


% --- distance from slab along flow path ---
if Vmax>=0
    gridD = gridSLC_vessel - gridSLC_vessel(:,:,1);
else
    gridD = gridSLC_vessel - gridSLC_vessel(:,:,end);
end
gridD = gridD./pVessel.n_hat(3);

% --- number of RF pulses ---
gridN = inf(size(vMap));
gridN(lumenMask) = gridD(lumenMask)/10./vMap(lumenMask)./pMri.TR;


% --- signal magnitude ---
mMap = zeros(size(gridN));
% -blood-
[Mz , pMri] = getMz_n(pMri, pMri.relax.blood, gridN(lumenMask));
[mMap(lumenMask), pMri] = getMxy_ss(Mz, pMri, pMri.relax.blood       );
% -tissue-
if ~isempty(pVessel.S.surround)
    mMap(tissueMask) = pVessel.S.surround;
    % S_surr = pVessel.S.surround;   % user-specified (e.g. S_tissue from fit)
else
    [Mz, pMri] = getMz_n(pMri, pMri.relax.blood, gridN(tissueMask));
    [mMap(tissueMask), pMri] = getMxy_ss(Mz, pMri, pMri.relax.GM);
end
% -wall-
if ~isempty(pVessel.S.wall)
    mMap(wallMask) = pVessel.S.wall;
else
    dbstack; error('code that')
end








% tmp = gridVoxIdxRaw;
% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorPE,squeeze(tmp(round(end/2),:,:))); axis image; colorbar
% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorFE,squeeze(tmp(:,round(end/2),:))); axis image; colorbar
% imagesc(pSim.spinGrid.coorPE,pSim.spinGrid.coorFE ,squeeze(tmp(:,:,round(end/2)))); axis image; colorbar







% % Handle 'cylinder3D' profile: fully 3D geometry with per-isochromat inflow signal.
% % Returns early so the 2D-only machinery below is untouched.
% if ischar(pVessel.geometry) && strcmp(pVessel.geometry, 'infiniteCylinder')
%     [magMap, vMap, pVessel, pSim, pMri] = simCylinder3D(pVessel, pSim, pMri);
%     return;
% end

% % ---- original simVesselSpins code (unchanged below this line) ----

% % Define radial coordinates (relative to vessel center)
% [gridFE, gridPE] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);  % dim1=FE, dim2=PE
% rGrid = sqrt((gridFE - pVessel.posFE).^2 + (gridPE - pVessel.posPE).^2);

% % Define compartment masks
% % Elliptical vessel if pVessel.AR and pVessel.alpha are provided; circular otherwise.
% R_lumen = pVessel.ID / 2;  % semi-major axis [mm]
% if isfield(pVessel,'AR') && ~isempty(pVessel.AR) && pVessel.AR ~= 1
%     alpha_v = pVessel.alpha;
%     AR_v    = pVessel.AR;
%     dFE = gridFE - pVessel.posFE;
%     dPE = gridPE - pVessel.posPE;
%     u   =  dPE.*cos(alpha_v) + dFE.*sin(alpha_v);   % along major axis
%     w   = -dPE.*sin(alpha_v) + dFE.*cos(alpha_v);   % along minor axis
%     ellipseDist = sqrt(u.^2 + (AR_v .* w).^2);      % = R_lumen at ellipse boundary
%     pVessel.mask.lumen    = ellipseDist <= R_lumen;
%     pVessel.mask.wall     = ellipseDist >  R_lumen & ellipseDist <= R_lumen + pVessel.WT;
%     pVessel.mask.surround = ~pVessel.mask.lumen & ~pVessel.mask.wall;
% else
%     pVessel.mask.lumen    = rGrid <= R_lumen;
%     pVessel.mask.wall     = rGrid >  R_lumen & rGrid <= R_lumen + pVessel.WT;
%     pVessel.mask.surround = rGrid >  R_lumen + pVessel.WT;
% end
% assert(all(pVessel.mask.lumen(:) + pVessel.mask.wall(:) + pVessel.mask.surround(:) == 1), ...
%     'pVessel.mask: lumen, wall, surround must be non-overlapping and cover all spins');

% % Define spin velocity map
% if ischar(pVessel.profile)
%     vMap = getVelMap(rGrid, pVessel.ID, pVessel.profile, pVessel.PD); % [cm/s]
%     if ~isempty(pVessel.vMax) && isempty(pVessel.vMean)
%         vMap = scale2maxVel(vMap, pVessel.vMax); % to the desired maximum velocity
%     elseif ~isempty(pVessel.vMean) && isempty(pVessel.vMax)
%         vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
%     elseif strcmp(pVessel.profile,'parabolic1') && ~isempty(pVessel.vMax) && ~isempty(pVessel.vMean) && pVessel.vMax/2==pVessel.vMean
%         vMap = scale2meanVel(vMap, pVessel.vMean, pVessel.mask.lumen); % to the desired mean velocity
%     else
%         error('Either pVessel.vMax or pVessel.vMean must be specified');
%     end
% elseif isnumeric(pVessel.profile)
%     vMap = zeros(size(rGrid));
%     vMap(:) = pVessel.profile;
% else
%     dbstack; error('Invalid vessel profile');
% end


% % MR signal magnitude
% % vessel lumen signal (flowing)
% if isempty(pVessel.S.lumen)
%     switch pVessel.profile
%         case 'plug'
%             [Mz_vMean,pMri] = getMz_ss(          pMri,pMri.relax.blood,pVessel.vMean);
%             [Mxy_vMax,pMri] = getMxy_ss(Mz_vMean,pMri,pMri.relax.blood              );
%             pVessel.S.lumen = Mxy_vMax;
%         case {'parabolic','parabolic1'}
%             [Mz ,pMri] = getMz_ss(    pMri,pMri.relax.blood,vMap(pVessel.mask.lumen));
%             [Mxy,pMri] = getMxy_ss(Mz,pMri,pMri.relax.blood                         );
%             pVessel.S.lumen = Mxy;
%         otherwise
%             dbstack; error('Invalid vessel profile');
%     end
% end
% % vessel surround (static)
% if isempty(pVessel.S.surround)
%     pVessel.S.surround = getMxy_ss(getMz_ss(pMri,pMri.relax.GM),pMri,pMri.relax.GM);
% end
% % map signal magnitude (divide by nSpinPerVox so summing spins within voxels gives S, the hypothetical measured signal if the voxel was single-compartment)
% magMap = zeros(size(rGrid));
% magMap(pVessel.mask.lumen)    = pVessel.S.lumen    ./pSim.isochromatPerVox;
% magMap(pVessel.mask.wall)     = pVessel.S.wall     ./pSim.nSpinPerVox;
% magMap(pVessel.mask.surround) = pVessel.S.surround ./pSim.nSpinPerVox;


% % Precompute montecarlo tessalation
% if pSim.monteCarloN > 0 && (~isfield(pSim,'monteCarloShiftFE') || ~isfield(pSim,'monteCarloShiftPE') || isempty(pSim.monteCarloShiftFE) || isempty(pSim.monteCarloShiftPE))
%     nSpinFE = pSim.spinGrid.matFE / pSim.voxGrid.matFE;  % spins per voxel in FE
%     shiftFE = (1:nSpinFE)-nSpinFE/2-0.5;
%     nSpinPE = pSim.spinGrid.matPE / pSim.voxGrid.matPE;  % spins per voxel in PE
%     shiftPE = (1:nSpinPE)-nSpinPE/2-0.5;
%     % find all possible combination of FE and PE shifts
%     [idx1, idx2] = ndgrid(1:length(shiftFE), 1:length(shiftPE));
%     idx = [idx1(:), idx2(:)];
%     % remove the no-shift combination since it is always done before
%     idx(all(idx==[find(shiftFE==0) find(shiftPE==0)],2),:) = [];
%     % shuffle
%     if pSim.monteCarloN~=inf
%         idx = idx(randperm(size(idx,1),pSim.monteCarloN),:);
%     end
%     pSim.monteCarloShiftFE = shiftFE(idx(:,1));
%     pSim.monteCarloShiftPE = shiftPE(idx(:,2));
% end


% % =========================================================================
% % Local functions
% % =========================================================================

% function [magMap, vMap, pVessel, pSim, pMri] = simCylinder3D(pVessel, pSim, pMri)
% % Full 3D cylinder geometry for the 'cylinder3D' profile.
% %
% % The vessel is an infinite cylinder defined by:
% %   n_hat  = pVessel.n_hat   [nx; ny; nz]  unit axis vector
% %   posFE  = pVessel.posFE   [mm]           FE center offset
% %   posPE  = pVessel.posPE   [mm]           PE center offset
% %   posSLC = pVessel.posSLC  [mm]           slice-direction offset (default 0)
% %   ID     = pVessel.ID      [mm]           lumen inner diameter
% %   WT     = pVessel.WT      [mm]           wall thickness (OD/2 = ID/2 + WT)
% %   Vmax   = pVessel.Vmax    [cm/s]         peak parabolic velocity
% %   A      = pVessel.A       [a.u.]         lumen signal amplitude scale
% %
% % Inflow model: instantaneous Mz after n_k pulses (per-isochromat),
% % computed via getMz_n + getMxy_ss. No effective-slice-thickness approximation.

% % n_hat  = pVessel.n_hat;   % [nx; ny; nz]
% % nx = n_hat(1);   ny = n_hat(2);   nz = n_hat(3);
% % % if isfield(pVessel,'posSLC');  cx_SLC = pVessel.posSLC;  else;  cx_SLC = 0;  end

% % % --- Z coordinate grid for the slice direction ---
% % % Sample Z at the same count as sub-spins per voxel in FE (e.g. 7 for a 7×7 grid).
% % % This is coarser than the FE/PE spin spacing but sufficient to resolve n_k(z).
% % % Slab thickness comes from pMri.sliceThickness; matSlice=1 assumed.
% % d_slab_zg  = pMri.sliceThickness;   % [mm]
% % nSpSLC     = round(numel(pSim.spinGrid.coorFE) / numel(pSim.voxGrid.coorFE));
% % nSpSLC     = max(1, nSpSLC);
% % if mod(nSpSLC, 2) == 0; nSpSLC = nSpSLC + 1; end   % make odd
% % dSLC       = d_slab_zg / nSpSLC;
% % coorZ      = ((-(nSpSLC-1)/2) : ((nSpSLC-1)/2)) .* dSLC;   % cell-centered [mm]

% % % Effective spin count per voxel including Z dimension
% % nSpinPerVox_3D = pSim.nSpinPerVox * nSpSLC;

% % --- 3D spin coordinate grids ---
% % Sizes: [nTotalFE, nTotalPE, nSpSLC]
% [gridFE, gridPE, gridSLC] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE, pSim.spinGrid.coorSLC);

% % --- 3D cylinder gridR ---
% fe_rel  = gridFE  - pVessel.pos(1);
% pe_rel  = gridPE  - pVessel.pos(2);
% slc_rel = gridSLC - pVessel.pos(3);
% % z_rel  = gridZ  - cx_SLC;
% dot_n  = fe_rel.*pVessel.n_hat(1) + pe_rel.*pVessel.n_hat(2) + slc_rel.*pVessel.n_hat(3);
% gridR = sqrt(max(0, fe_rel.^2 + pe_rel.^2 + slc_rel.^2 - dot_n.^2));
% % imagesc(pSim.spinGrid.coorPE,pSim.spinGrid.coorFE,gridR(:,:,round(end/2))); axis image

% % --- Compartment masks (3D) ---
% R_blood = pVessel.ID / 2;
% R_wall  = R_blood + pVessel.WT;
% lumen3  = gridR <= R_blood;
% wall3   = gridR >  R_blood & gridR <= R_wall;
% surr3   = ~lumen3 & ~wall3;

% % --- Parabolic velocity in cylinder coordinates ---
% switch pVessel.profile
%     case 'parabolic1'
%         Vmax  = pVessel.vMean*2;
%     otherwise
%         dbstack; error('double-check that')
% end
% A_scl = pSim.mriScale;
% v3    = zeros(size(gridR));
% v3(lumen3) = Vmax .* (1 - (gridR(lumen3) ./ R_blood).^2);

% % --- Per-isochromat number of RF pulses ---
% gridD = (slc_rel + pMri.sliceThickness/2) ./ pVessel.n_hat(3);   % [mm] path from slab entry
% gridN = ceil(gridD/10./v3./pMri.TR);

% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorPE,squeeze(gridD(round(end/2),:,:))); axis image
% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorFE,squeeze(gridD(:,round(end/2),:))); axis image
% imagesc(pSim.spinGrid.coorPE,pSim.spinGrid.coorFE,squeeze(gridD(:,:,round(end/2)))); axis image


% % gridN = ceil(max(0, gridD/10) ./ v3 ./ pMri.TR);
% % gridN(v3 <= 0) = Inf;   % stationary spins → steady-state saturation





% % d_slab = pMri.sliceThickness;   % [mm]  true slab thickness

% % % Path along cylinder axis from slab entry (Z = -d_slab/2 + cx_SLC) to current Z position.
% % % Assumes nz > 0 (blood flows in +z direction through the slab).

% % z_lumen    = slc_rel(lumen3);


% % d_trav     = (z_lumen + d_slab/2) ./ max(pVessel.n_hat(3), eps);   % [mm] path from slab entry
% % v_lumen    = v3(lumen3);                              % [cm/s]
% % % Convert d_trav mm → cm before dividing by v [cm/s] and TR [s].
% % n_k        = ceil(max(0, d_trav/10) ./ max(v_lumen, eps) ./ pMri.TR);
% % n_k(v_lumen < eps) = Inf;   % stationary spins → steady-state saturation

% % getMz_n: instantaneous Mz after exactly n_k pulses from fully-relaxed M0.
% [Mz_k, pMri] = getMz_n(pMri, pMri.relax.blood, gridN);
% % getMxy_ss: apply sin(FA) flip and T2* decay to get Mxy.
% [Mxy_k, pMri] = getMxy_ss(Mz_k, pMri, pMri.relax.blood);

% % --- Surround signal ---
% if ~isempty(pVessel.S.surround)
%     S_surr = pVessel.S.surround;   % user-specified (e.g. S_tissue from fit)
% else
%     % Auto-compute from GM relaxation (default for non-fitting use)
%     [Mz_gm, pMri]  = getMz_ss(pMri, pMri.relax.GM, 0);
%     [S_surr, pMri] = getMxy_ss(Mz_gm, pMri, pMri.relax.GM);
% end

% % --- Build 3D signal array and Z-sum to 2D ---
% % S3: [nTotalFE, nTotalPE, nSpSLC]
% S3 = zeros(size(gridR));
% S3(lumen3) = A_scl .* Mxy_k(lumen3);   % per-lumen-spin inflow signal
% S3(surr3)  = S_surr;           % scalar, broadcast over surround spins
% % wall3 contribution = 0 (already zero)

% % % Sum over Z dimension → [nTotalFE, nTotalPE]
% % S2 = sum(S3, 3);

% % --- 3D masks and radial map (full isochromat resolution, no projection) ---
% assert(all(lumen3(:) + wall3(:) + surr3(:) == 1), ...
%     'simCylinder3D: 3D masks must be non-overlapping and exhaustive');

% pVessel.mask.lumen    = lumen3;               % [nTotalFE, nTotalPE, nSpSLC]
% pVessel.mask.wall     = wall3;
% pVessel.mask.surround = surr3;
% pVessel.rGrid         = single(gridR);       % [nTotalFE, nTotalPE, nSpSLC] radial distance from axis [mm]
% pVessel.dGrid         = single(gridD);       
% pVessel.nGrid         = single(gridN);       
% pVessel.S.lumen       = single(S3(lumen3));   % per-isochromat inflow signal
% pVessel.S.wall        = 0;
% pVessel.S.surround    = S_surr;

% % --- 2D magMap: Z-sum of 3D signal, normalised by 3D spin count ---
% magMap = sum(S3, 3) ./ nSpinPerVox_3D;   % [nTotalFE, nTotalPE]

% % --- 2D velocity map (at Z=0, slab midplane) ---
% fe_2d = gridFE(:,:,1) - pVessel.posFE;
% pe_2d = gridPE(:,:,1) - pVessel.posPE;
% z_2d  = 0 - cx_SLC;   % scalar
% dot_n_2d  = fe_2d.*nx + pe_2d.*ny + z_2d.*nz;
% gridR_2d = sqrt(max(0, fe_2d.^2 + pe_2d.^2 + z_2d.^2 - dot_n_2d.^2));
% vMap  = Vmax .* max(0, 1 - (gridR_2d ./ R_blood).^2);
% vMap(gridR_2d > R_blood) = 0;
