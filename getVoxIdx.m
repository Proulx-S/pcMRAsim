function voxIdx = getVoxIdx(voxGrid, spinGrid)
% getVoxIdx  Distance-sorted spin-to-voxel index map.
%   voxIdx = getVoxIdx(voxGrid, spinGrid)
%
%   Returns a (spinGrid.matPE x spinGrid.matFE) matrix of voxel indices.
%   Index 0 = voxel whose center is closest to (0,0); increasing indices
%   for voxels further away. Ties broken by FE center coord, then PE.

[gridFE, gridPE] = ndgrid(spinGrid.coorFE, spinGrid.coorPE);  % dim1=FE, dim2=PE

iFE = floor(gridFE / voxGrid.dFE + voxGrid.matFE/2) + 1;
iPE = floor(gridPE / voxGrid.dPE + voxGrid.matPE/2) + 1;
iFE = max(1, min(voxGrid.matFE, iFE));
iPE = max(1, min(voxGrid.matPE, iPE));
gridVoxIdxRaw = sub2ind([voxGrid.matFE, voxGrid.matPE], iFE, iPE);

[voxFE_grid, voxPE_grid] = ndgrid(1:voxGrid.matFE, 1:voxGrid.matPE);
voxCtrFE = voxGrid.coorFE(voxFE_grid);  % [mm]
voxCtrPE = voxGrid.coorPE(voxPE_grid);  % [mm]
voxDist  = sqrt(voxCtrFE.^2 + voxCtrPE.^2);

[~, sortOrder] = sortrows([voxDist(:), voxCtrFE(:), voxCtrPE(:)]);
newIdxMap = zeros(voxGrid.matFE * voxGrid.matPE, 1);
newIdxMap(sortOrder) = 0 : voxGrid.matFE*voxGrid.matPE-1;
voxIdx = reshape(newIdxMap(gridVoxIdxRaw(:)), size(gridVoxIdxRaw));
