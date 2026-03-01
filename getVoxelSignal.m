function [I,voxIdx] = getVoxelSignal(spinMap,pSim,mask)

% Define voxel indices
if isstruct(pSim)
    nX     = pSim.FOVx./pSim.voxSizeX;
    nY     = pSim.FOVy./pSim.voxSizeY;
    [nRows, nCols] = size(spinMap);
    [xBin, ~] = discretize(1:nCols, round(linspace(1, nCols+1, nX+1)));
    [yBin, ~] = discretize(1:nRows, round(linspace(1, nRows+1, nY+1)));
    [xx, yy] = meshgrid(xBin, yBin);  % xx: x-bin index per column, yy: y-bin index per row
    voxIdx = (xx - 1) * nY + yy;
else
    voxIdx = pSim;
    nX     = length(unique(voxIdx(1,:)));
    nY     = length(unique(voxIdx(:,1)));
end

% Get voxel signal
I = zeros(nY,nX);
for iVoxel = 1:(nX*nY)
    I(iVoxel) = sum(spinMap(voxIdx==iVoxel & mask));
end