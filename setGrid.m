% function [voxGrid, spinGrid, nSpinPerVox] = setGrid(fovFE, fovPE, matFE, matPE, nSpin, mode)
function [voxGrid, spinGrid, isochromatPerVox] = setGrid(fov, mat, isochromatPerVoxPerDirMin, mode)
% INPUTS
%  fovFE, fovPE              : field of view [mm]
%  matFE, matPE              : number of voxels
%  isochromatPerVoxPerDirMin : minimum number of isochromat in each direction in a voxel [n]
%  mode                      : 'pseudoVoxel' (default) | 'centerVox'
% OUTPUTS
%  voxGrid     — voxel grid (one entry per voxel)
%    .fovFE, .fovPE   [mm]     total FOV
%    .matFE, .matPE   [—]      number of voxels
%    .dFE,   .dPE     [mm]     voxel size
%    .coorFE, .coorPE [mm]     1D voxel-center coordinate vectors
%  spinGrid    — spin grid (one entry per spin total)
%    .fovFE, .fovPE   [mm]     total FOV  (same as voxGrid)
%    .matFE, .matPE   [—]      total spins in each direction (nSpinPerVox * voxGrid.mat)
%    .dFE,   .dPE     [mm]     spin spacing
%    .coorFE, .coorPE [mm]     1D spin coordinate vectors (full grid)
%  isochromatPerVox   [n]      isochromat per voxel
%
% Recoverable quantities (not stored):
%   nSpin total                = spinGrid.matFE * spinGrid.matPE
%   voxel boundaries FE        = voxGrid.coorFE +/- voxGrid.dFE/2
%   spin-to-voxel map          = getVoxIdx(voxGrid, spinGrid)
%   2D spin coord grids        : [spinGridFE,spinGridPE] = ndgrid(spinGrid.coorFE, spinGrid.coorPE)  [dim1=FE, dim2=PE]
%   2D spin radial coord grids : voxGridR = sqrt(voxGridFE.^2+voxGridPE.^2);  [dim1=FE, dim2=PE]

if ~exist('mode','var') || isempty(mode); mode = 'pseudoVoxel'; end

if strcmp(mode,'centerVox')
    assert(mod(matFE,2)==1 && mod(matPE,2)==1 && mod(matSLC,2)==1, 'centerVox mode requires odd matFE and matPE');
end

% fov    = [fovFE fovPE fovSLC];
% mat = [matFE matPE matSLC];
voxSz  = fov./mat; % [mm]

% % Voxel size [mm]
% voxSzFE   = fovFE  / matFE;    % [mm]
% voxSzPE   = fovPE  / matPE;    % [mm]
% voxSzSLC  = fovSLC / matSLC;   % [mm]
% % fovAspect = fovFE  / fovPE;    % [—]

% Odd number of spins per voxel in each direction, preserving FOV aspect ratio
% (equal physical spin spacing: dFE ≈ dPE ≈ sqrt(fovFE*fovPE/nSpin) [mm])
[~,b] = sort(voxSz,'ascend');
isochromatPerVoxPerDirMin = isochromatPerVoxPerDirMin + mod(isochromatPerVoxPerDirMin+1, 2);  % make odd
d = voxSz(b(1))/isochromatPerVoxPerDirMin;
nIsochromat = nan(size(voxSz));
nIsochromat(b==1) = voxSz(b==1)./d;
nIsochromat(b~=1) = round(voxSz(b~=1)./d);
nIsochromat = nIsochromat + mod(nIsochromat+1, 2);  % make odd
d           = voxSz./nIsochromat;    % [mm]
nTotal      = nIsochromat.*mat ;


% nSpinFE   = ceil(voxSzFE /d);
% nSpinPE   = ceil(voxSzPE /d);
% nSpinSLC  = ceil(voxSzSLC/d);

% d = nthroot(voxSzFE*voxSzPE*voxSzSLC/isochromatPerVoxPerDirMin,3);
% nSpinFE  = round(voxSzFE/d);
% nSpinFE  = nSpinFE + mod(nSpinFE+1, 2);  % make odd
% nSpinPE  = round(voxSzPE/d);
% nSpinPE  = nSpinPE + mod(nSpinPE+1, 2);  % make odd
% nSpinSLC = round(voxSzSLC/d);
% nSpinSLC = nSpinSLC + mod(nSpinSLC+1, 2);  % make odd
% % nSpinFE = round(sqrt(nSpin * fovAspect) / matFE);
% % nSpinFE = nSpinFE + mod(nSpinFE+1, 2);  % make odd
% % nSpinPE = round(sqrt(nSpin / fovAspect) / matPE);
% % nSpinPE = nSpinPE + mod(nSpinPE+1, 2);  % make odd

% % Spin spacing [mm]
% dFE  = voxSzFE  / nSpinFE;    % [mm]
% dPE  = voxSzPE  / nSpinPE;    % [mm]
% dSLC = voxSzSLC / nSpinSLC;    % [mm]

% % Total spins per direction
% nTotalFE  = nSpinFE  * matFE ;
% nTotalPE  = nSpinPE  * matPE ;
% nTotalSLC = nSpinSLC * matSLC;

% Voxel grid
voxGrid.fovFE  = fov(1);
voxGrid.fovPE  = fov(2);
voxGrid.fovSLC = fov(3);
voxGrid.matFE  = mat(1);
voxGrid.matPE  = mat(2);
voxGrid.matSLC = mat(3);
voxGrid.dFE    = voxSz(1) ;    % [mm]
voxGrid.dPE    = voxSz(2) ;    % [mm]
voxGrid.dSLC   = voxSz(3) ;    % [mm]
voxGrid.coorFE  = (-(mat(1)-1)/2 : (mat(1)-1)/2) * voxSz(1);  % [mm]  1 x matFE
voxGrid.coorPE  = (-(mat(2)-1)/2 : (mat(2)-1)/2) * voxSz(2);  % [mm]  1 x matPE
voxGrid.coorSLC = (-(mat(3)-1)/2 : (mat(3)-1)/2) * voxSz(3);  % [mm]  1 x matSLC

% Spin grid  (fovFE = matFE * dFE holds for both grids by construction)
spinGrid.fovFE  = fov(1);
spinGrid.fovPE  = fov(2);
spinGrid.fovSLC = fov(3);
spinGrid.matFE  = nTotal(1);
spinGrid.matPE  = nTotal(2);
spinGrid.matSLC = nTotal(3);
spinGrid.dFE    = d(1);        % [mm]
spinGrid.dPE    = d(2);        % [mm]
spinGrid.dSLC   = d(3);        % [mm]
spinGrid.coorFE  = (-(nTotal(1)-1)/2 : (nTotal(1)-1)/2) * d(1);  % [mm]  1 x nTotalFE
spinGrid.coorPE  = (-(nTotal(2)-1)/2 : (nTotal(2)-1)/2) * d(2);  % [mm]  1 x nTotalPE
spinGrid.coorSLC = (-(nTotal(3)-1)/2 : (nTotal(3)-1)/2) * d(3);  % [mm]  1 x nTotalSLC

isochromatPerVox = prod(nTotal./mat);  % spins per voxel
