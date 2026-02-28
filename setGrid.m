function [xGrid, yGrid, rGrid, dx, dy, nSpin, spinDensity] = setGrid(FOVx, FOVy, nSpin)
% FOVx :  Field of view in x-direction [mm]
% FOVy :  Field of view in y-direction [mm]
% dx   :  grid spacing in x-direction [mm]
% dy   :  grid spacing in y-direction [mm]
% nSpin:  If length(nSpin)=1 -> desired approximate total number of spins in the grid
% nSpin:  If length(nSpin)=2 -> desired number of rows x columns of the spin grid
% nSpin:  as an output, it is of the length(nSpin)=1 format
% xGrid:  x-coordinates of the grid [mm]
% yGrid:  y-coordinates of the grid [mm]
% spinDensity:  Number of spins per mm^2

% Handle input nSpin: can be scalar (total spins) or vector (rows, cols)
if length(nSpin) == 1
    % Approximate total number of spins desired
    % Calculate grid dimensions to achieve roughly nSpin total spins
    targetDensity = nSpin / (FOVx * FOVy); % spins/mm^2
    xN = round(FOVx * sqrt(targetDensity)); % number of spins in x-direction
    yN = round(FOVy * sqrt(targetDensity)); % number of spins in y-direction
elseif length(nSpin) == 2
    % Exact number of rows and columns specified
    yN = nSpin(1); % rows
    xN = nSpin(2); % columns
else
    error('nSpin must be a scalar or a 2-element vector');
end

% Create grid coordinates
x = linspace(-FOVx/2 + FOVx/(2*xN), FOVx/2 - FOVx/(2*xN), xN);
y = linspace(-FOVy/2 + FOVy/(2*yN), FOVy/2 - FOVy/(2*yN), yN);
[xGrid, yGrid] = meshgrid(x, y);

% Calculate actual spacing
dx = mean(diff(xGrid(1,:)));
dy = mean(diff(yGrid(:,1)));

% Calculate actual spin density
spinDensity = 1 / (dx * dy); % spins/mm^2

% Return total number of spins
nSpin = xN * yN;

% Calculate radial distance from the center of the grid [mm]
rGrid = sqrt(xGrid.^2 + yGrid.^2); % radial distance from the center of the grid [mm]