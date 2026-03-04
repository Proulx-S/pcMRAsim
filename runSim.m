function res = runSim(pVessel, pVenc, pSim, verbose, anaFlag, vFix)
if ~exist('verbose','var') || isempty(verbose); verbose = true; end
if ~exist('anaFlag','var')  || isempty(anaFlag); anaFlag = 'inflowOnSpinVelocity'; end
if ~exist('vFix','var')     || isempty(vFix);     vFix = []; end
    % anaFlag = 'inflowFixedAtMax';
    % anaFlag = 'inflowFixedAtVelocity';
    % anaFlag = 'inflowOnMeanVelocity';
    % anaFlag = 'inflowOnSpinVelocity';

% Define simulation grid
FOVx = pSim.FOVx;
FOVy = pSim.FOVy;
nVOXx = ceil(FOVx./pSim.voxSizeX);
nVOXy = ceil(FOVy./pSim.voxSizeY);
[xGrid, yGrid,   ~, pSim.dx, pSim.dy, pSim.nSpin, pSim.spinDensity] = setGrid(FOVx, FOVy, nVOXx, nVOXy, pSim.nSpin);


% Define velocity encoding
pVenc.mVenc; % velocity encoding bandwidth [cm/s]
pVenc.dVenc; % velocity encoding resolution [cm/s]
if pVenc.dVenc==inf
    venc = pVenc.mVenc;
    M1   = vencToM1(venc);
else
    M1   = linspace(0, vencToM1(pVenc.dVenc), round(vencToM1(pVenc.dVenc)/vencToM1(pVenc.mVenc))+1);
    venc = M1toVenc(M1);
end
pVenc.vencList = venc;


% Get maps and voxel signals for vessel at (x0,y0)
if verbose; disp('Vessel at (x0,y0). Simulating...'); end
[magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, pSim, [], [], anaFlag, vFix);
nX = pSim.FOVx./pSim.voxSizeX;
nY = pSim.FOVy./pSim.voxSizeY;
If = zeros(nY,nX,length(venc));
Is = zeros(nY,nX,length(venc));
for iVenc = 1:length(venc)
    % Apply velocity encoding to spins
    spinMap  = magMap.*exp(1i*vel2phase(vMap, venc(iVenc)));
    % Get voxel signal for each compartment
    [If(:,:,iVenc),voxIdx] = getVoxelSignal(spinMap,pSim  ,mask.lumen   );
    [Is(:,:,iVenc),     ~] = getVoxelSignal(spinMap,voxIdx,mask.surround);
end
res.If = If;
res.Is = Is;
res.pVessel = pVessel;
res.pVenc   = pVenc;
res.pSim    = pSim;
res.mask    = mask;
res.vMap    = vMap;
res.magMap  = magMap;
if verbose; disp('Vessel at (x0,y0). Done.'); end




if pSim.voxRndOffset
    % Randomize vessel position relative to center voxel of FOV
    xSpinRange = (pSim.voxSizeX./pSim.dx-1)/2 .*[-1 1];
    ySpinRange = (pSim.voxSizeY./pSim.dy-1)/2 .*[-1 1];
    x0 = randi(xSpinRange,[1 pSim.voxRndOffset]).*pSim.dx;
    y0 = randi(ySpinRange,[1 pSim.voxRndOffset]).*pSim.dy;

    % Get maps and voxel signals at randomized vessel positions
    if verbose; disp('Vessel at random positions. Simulating...'); end
    If = zeros(nY,nX,length(venc),length(x0));
    Is = zeros(nY,nX,length(venc),length(x0));
    res.mask.lumenRnd = zeros(size(xGrid));
    for iPos = 1:length(x0)
        % Define vessel maps
        [magMap,vMap,mask,pVessel] = simVesselSpins(xGrid, yGrid, pVessel, pSim, x0(iPos), y0(iPos), anaFlag, vFix);
        res.mask.lumenRnd = res.mask.lumenRnd+mask.lumen;

        for iVenc = 1:length(venc)
            % Apply velocity encoding to spins
            spinMap  = magMap.*exp(1i*vel2phase(vMap, venc(iVenc)));
            % Get voxel signal for each compartment
            [If(:,:,iVenc,iPos),voxIdx] = getVoxelSignal(spinMap,pSim  ,mask.lumen   );
            [Is(:,:,iVenc,iPos),     ~] = getVoxelSignal(spinMap,voxIdx,mask.surround);
        end
    end

    res.rndOffset.If = If;
    res.rndOffset.Is = Is;
    if verbose; disp('Vessel at random positions. Done.'); end
end
