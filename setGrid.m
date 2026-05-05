% function [voxGrid, spinGrid, nSpinPerVox] = setGrid(fovFE, fovPE, matFE, matPE, nSpin, mode)
function [voxGrid, spinGrid, nSpinPerVox] = setGrid(fovFE, fovPE, matFE, matPE, nSpin, mode)
% INPUTS
%  fovFE, fovPE : field of view [mm]
%  matFE, matPE : number of voxels
%  nSpin        : target total spin count [n]
%  mode         : 'pseudoVoxel' (default) | 'centerVox'
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
%  nSpinPerVox [n]             spins per voxel
%
% Recoverable quantities (not stored):
%   nSpin total                = spinGrid.matFE * spinGrid.matPE
%   voxel boundaries FE        = voxGrid.coorFE +/- voxGrid.dFE/2
%   spin-to-voxel map          = getVoxIdx(voxGrid, spinGrid)
%   2D spin coord grids        : [spinGridFE,spinGridPE] = ndgrid(spinGrid.coorFE, spinGrid.coorPE)  [dim1=FE, dim2=PE]
%   2D spin radial coord grids : voxGridR = sqrt(voxGridFE.^2+voxGridPE.^2);  [dim1=FE, dim2=PE]

if ~exist('mode','var') || isempty(mode); mode = 'pseudoVoxel'; end

if strcmp(mode,'centerVox')
    assert(mod(matFE,2)==1 && mod(matPE,2)==1, 'centerVox mode requires odd matFE and matPE');
end

% Voxel size [mm]
voxSzFE   = fovFE / matFE;    % [mm]
voxSzPE   = fovPE / matPE;    % [mm]
fovAspect = fovFE / fovPE;    % [—]

% Odd number of spins per voxel in each direction, preserving FOV aspect ratio
% (equal physical spin spacing: dFE ≈ dPE ≈ sqrt(fovFE*fovPE/nSpin) [mm])
nSpinFE = round(sqrt(nSpin * fovAspect) / matFE);
nSpinFE = nSpinFE + mod(nSpinFE+1, 2);  % make odd
nSpinPE = round(sqrt(nSpin / fovAspect) / matPE);
nSpinPE = nSpinPE + mod(nSpinPE+1, 2);  % make odd

% Spin spacing [mm]
dFE = voxSzFE / nSpinFE;    % [mm]
dPE = voxSzPE / nSpinPE;    % [mm]

% Total spins per direction
nTotalFE = nSpinFE * matFE;
nTotalPE = nSpinPE * matPE;

% Voxel grid
voxGrid.fovFE  = fovFE;
voxGrid.fovPE  = fovPE;
voxGrid.matFE  = matFE;
voxGrid.matPE  = matPE;
voxGrid.dFE    = voxSzFE;    % [mm]
voxGrid.dPE    = voxSzPE;    % [mm]
voxGrid.coorFE = (-(matFE-1)/2 : (matFE-1)/2) * voxSzFE;  % [mm]  1 x matFE
voxGrid.coorPE = (-(matPE-1)/2 : (matPE-1)/2) * voxSzPE;  % [mm]  1 x matPE

% Spin grid  (fovFE = matFE * dFE holds for both grids by construction)
spinGrid.fovFE  = fovFE;
spinGrid.fovPE  = fovPE;
spinGrid.matFE  = nTotalFE;
spinGrid.matPE  = nTotalPE;
spinGrid.dFE    = dFE;        % [mm]
spinGrid.dPE    = dPE;        % [mm]
spinGrid.coorFE = (-(nTotalFE-1)/2 : (nTotalFE-1)/2) * dFE;  % [mm]  1 x nTotalFE
spinGrid.coorPE = (-(nTotalPE-1)/2 : (nTotalPE-1)/2) * dPE;  % [mm]  1 x nTotalPE

nSpinPerVox = nSpinFE * nSpinPE;  % spins per voxel
