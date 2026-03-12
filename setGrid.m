function [gridFE, gridPE, gridVoxIdx, dFE, dPE, nSpin] = setGrid(voxSzFE, voxSzPE, matSzFE, matSzPE, nSpin)
% INPUTS
%  voxSzFE:  voxel size in FE direction [mm]
%  voxSzPE:  voxel size in PE direction [mm]
%  matSzFE:  matrix size in FE direction [voxels]
%  matSzPE:  matrix size in PE direction [voxels]
%  nSpin  :  desired approximate number of spins per voxel [n]
% OUTPUTS
%  gridFE    :  FE grid coordinates [mm]
%  gridPE    :  PE grid coordinates [mm]
%  gridVoxIdx:  voxel indices of spins with idx=0 for center voxel
%  dFE       :  spin spacing in FE direction [mm]
%  dPE       :  spin spacing in PE direction [mm]
%  nSpin     :  actual number of spins in a voxel [n]


% odd number of spins in each directions in a voxel
nSpinFE = round(sqrt(nSpin*voxSzFE/voxSzPE));
nSpinPE = round(sqrt(nSpin*voxSzPE/voxSzFE));
nSpinFE = nSpinFE + mod(nSpinFE+1,2);
nSpinPE = nSpinPE + mod(nSpinPE+1,2);
% spin spacing in each directions
dFE     = voxSzFE/nSpinFE;
dPE     = voxSzPE/nSpinPE;
% number of spins in a voxel
nSpin   = nSpinFE * nSpinPE;

% spin cartesian coordinates relative to center of center voxel
gridFE = linspace(-voxSzFE*matSzFE/2+dFE/2, voxSzFE*matSzFE/2-dFE/2, nSpinFE);
gridPE = linspace(-voxSzPE*matSzPE/2+dPE/2, voxSzPE*matSzPE/2-dPE/2, nSpinPE);
[gridFE, gridPE] = meshgrid(gridFE, gridPE);

% voxel indices of spins: matrix same size as gridR, one index per voxel
iFE = round(gridFE/voxSzFE + (matSzFE+1)/2);
iPE = round(gridPE/voxSzPE + (matSzPE+1)/2);
iFE = max(1, min(matSzFE, iFE));
iPE = max(1, min(matSzPE, iPE));
gridVoxIdx = sub2ind([matSzFE, matSzPE], iFE, iPE);
% make center voxel index 0
idx0 = gridVoxIdx(round(end/2),round(end/2));
gridVoxIdx(gridVoxIdx==idx0) = 0;
gridVoxIdx(gridVoxIdx> idx0) = gridVoxIdx(gridVoxIdx> idx0)-1;
