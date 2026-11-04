function voxIdx = getVoxIdx(voxGrid, spinGrid, gridFE, gridPE, gridSLC, gridR)
% getVoxIdx  Distance-sorted spin-to-voxel index map.
%   voxIdx = getVoxIdx(voxGrid, spinGrid)
%
%   Returns a (spinGrid.matPE x spinGrid.matFE) matrix of voxel indices.
%   Index 0 = voxel whose center is closest to (0,0); increasing indices
%   for voxels further away. Ties broken by FE center coord, then PE.

if ~exist('gridFE') || isempty(gridFE) || ~exist('gridPE') || isempty(gridPE) || ~exist('gridSLC') || isempty(gridSLC)
    dbstack; error('double-check that');
    [gridFE, gridPE] = ndgrid(spinGrid.coorFE, spinGrid.coorPE);  % dim1=FE, dim2=PE
end

iFE  = floor(gridFE  / voxGrid.dFE  + voxGrid.matFE /2) + 1;
iPE  = floor(gridPE  / voxGrid.dPE  + voxGrid.matPE /2) + 1;
iSLC = floor(gridSLC / voxGrid.dSLC + voxGrid.matSLC/2) + 1;
iFE  = max(1, min(voxGrid.matFE , iFE ));
iPE  = max(1, min(voxGrid.matPE , iPE ));
iSLC = max(1, min(voxGrid.matSLC, iSLC));
voxIdx = uint16(sub2ind([voxGrid.matFE, voxGrid.matPE, voxGrid.matSLC], iFE, iPE, iSLC));

% % Sort 
% [~, sortOrder] = sortrows([gridR(:), gridFE(:), gridPE(:), gridSLC(:)]);
% newIdxMap = zeros(voxGrid.matFE * voxGrid.matPE * voxGrid.matSLC, 1);
% newIdxMap(sortOrder) = 0 : voxGrid.matFE*voxGrid.matPE*voxGrid.matSLC-1;
% voxIdx = reshape(newIdxMap(gridVoxIdxRaw(:)), size(gridVoxIdxRaw));



% tmp = gridVoxIdxRaw;
% imagesc(spinGrid.coorSLC,spinGrid.coorPE,squeeze(tmp(round(end/2),:,:))); axis image; colorbar
% imagesc(spinGrid.coorSLC,spinGrid.coorFE,squeeze(tmp(:,round(end/2),:))); axis image; colorbar
% imagesc(spinGrid.coorPE,spinGrid.coorFE ,squeeze(tmp(:,:,round(end/2)))); axis image; colorbar


% [voxFE_grid, voxPE_grid, voxSLC_grid] = ndgrid(1:voxGrid.matFE, 1:voxGrid.matPE, 1:voxGrid.matSLC);
% voxCtrFE  = voxGrid.coorFE(voxFE_grid);  % [mm]
% voxCtrPE  = voxGrid.coorPE(voxPE_grid);  % [mm]
% voxCtrSLC = voxGrid.coorPE(voxSLC_grid);  % [mm]
% voxDist = sqrt(voxCtrFE.^2 + voxCtrPE.^2 + voxCtrSLC.^2);

% [~, sortOrder] = sortrows([voxDist(:), voxCtrFE(:), voxCtrPE(:)]);
% newIdxMap = zeros(voxGrid.matFE * voxGrid.matPE, 1);
% newIdxMap(sortOrder) = 0 : voxGrid.matFE*voxGrid.matPE-1;
% voxIdx = reshape(newIdxMap(gridVoxIdxRaw(:)), size(gridVoxIdxRaw));
