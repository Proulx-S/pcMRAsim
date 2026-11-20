% function [gridFE, gridPE, gridVoxIdx, dFE, dPE, nSpin] = setGrid(voxSzFE, voxSzPE, matFE, matPE, nSpin)
function [gridFE, gridPE, gridVoxIdx, dFE, dPE, nSpin] = setGrid(fovFE, fovPE, matFE, matPE, nSpin)
% INPUTS
%  fovFE:  field of view in FE direction [mm]
%  fovPE:  field of view in PE direction [mm]
%  matFE:  matrix size in FE direction [voxels]
%  matPE:  matrix size in PE direction [voxels]
%  nSpin  :  desired approximate number of spins per voxel [n]
% OUTPUTS
%  gridFE    :  FE grid coordinates [mm]
%  gridPE    :  PE grid coordinates [mm]
%  gridVoxIdx:  voxel indices of spins with idx=0 for center voxel
%  dFE       :  spin spacing in FE direction [mm]
%  dPE       :  spin spacing in PE direction [mm]
%  nSpin     :  actual number of spins in a voxel [n]


% odd number of spins in each directions in a voxel
voxSzFE = fovFE/matFE;
voxSzPE = fovPE/matPE;
voxAspectRatio = voxSzFE/voxSzPE;
nSpinFE = round(sqrt(nSpin*voxAspectRatio));
nSpinPE = round(sqrt(nSpin*voxAspectRatio));
nSpinFE = nSpinFE + mod(nSpinFE+1,2);
nSpinPE = nSpinPE + mod(nSpinPE+1,2);
% spin spacing in each directions
dFE     = voxSzFE/nSpinFE;
dPE     = voxSzPE/nSpinPE;
% number of spins in a voxel
nSpin   = nSpinFE * nSpinPE;

% spin cartesian coordinates relative to center of center voxel
gridFE = linspace(-fovFE/2+dFE/2, fovFE/2-dFE/2, nSpinFE);
gridPE = linspace(-fovPE/2+dPE/2, fovPE/2-dPE/2, nSpinPE);
[gridFE, gridPE] = meshgrid(gridFE, gridPE);

% voxel indices of spins: matrix same size as gridR, one index per voxel
iFE = round(gridFE/voxSzFE + (matFE+1)/2);
iPE = round(gridPE/voxSzPE + (matPE+1)/2);
iFE = max(1, min(matFE, iFE));
iPE = max(1, min(matPE, iPE));
gridVoxIdx = sub2ind([matFE, matPE], iFE, iPE);
% make center voxel index 0
idx0 = gridVoxIdx(round(end/2),round(end/2));
gridVoxIdx(gridVoxIdx==idx0) = 0;
gridVoxIdx(gridVoxIdx> idx0) = gridVoxIdx(gridVoxIdx> idx0)-1;
