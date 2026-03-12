function res = runSim(pVessel, pVenc, pSim, verbose)
if ~exist('verbose','var') || isempty(verbose); verbose = true; end

if nargin==0
    %% Output default parameters

    pMri.sliceThickness     = 1;      % [mm]   | default 1 mm
    pMri.TR                 = 0.05;   % [s]    | default 0.05 s
    pMri.FA                 = 40;     % [deg]  | default 40 deg
    pMri.venc.vencRes       = 2;      % [cm/s] | default 2 cm/s
    pMri.venc.vencMax       = 50;     % [cm/s] | default 50 cm/s
    pMri.venc.vencList      = M1toVenc( linspace(0, vencToM1(pMri.venc.vencRes), round(vencToM1(pMri.venc.vencRes)/vencToM1(pMri.venc.vencMax))+1) );
    % pMri.venc.vencList      = [inf 8]; % [cm/s] | default [inf 8] or a list built from vencRes and vencMax
    % pMri.venc.vencRes       = [];
    % pMri.venc.vencMax       = [];
    pMri.relax.blood.T1     = 2.58;   % [s]    | default human in vivo at 7T
    % Blood T1. Human at 7T (Rane & Gore, Magn Reson Imaging 31(3):477–479, 2013, doi:10.1016/j.mri.2012.08.008):
    %   arterial 2.29±0.10 s, venous 2.07±0.12 s in vitro (37°C); venous sagittal sinus in vivo 2.45±0.11 s.
    %   arterial in vivo 2.45/2.07*2.29 = 2.71 s
    %   mid arterio-venous in vivo (2.71+2.45)/2 = 2.58 s
    pMri.relax.blood.T2star = 0.01;   % [s]    | default human in vivo at 7T
    % Blood T2*. Human at 7T: venous blood ~7.4 ms (SWI/venography at 7T); arterial longer (higher oxygenation).
    %   Blood T2* is strongly oxygenation-dependent; R2* increases with deoxyhemoglobin (e.g. Qin & van Zijl, MRM 24868, 2009).
    pMri.relax.GM.T1        = 1.939;  % [s]    | default human in vivo at 7T
    % Gray matter cortical T1. Human at 7T (Waddell et al., MAGMA 21:121–130, 2008, doi:10.1007/s10334-008-0104-8):
    %   cortical gray matter 1939±149 ms (~1.94 s), white matter 1126±97 ms (MPRAGE, 4 subjects).
    pMri.relax.GM.T2star    = 0.0329; % [s]    | default human in vivo at 7T
    % Gray matter cortical T2*. Human at 7T (Peters et al., Magn Reson Imaging 25:748–753, 2007, doi:10.1016/j.mri.2007.02.014):
    %   cortical gray matter 32.9±2.3 ms, white matter 27.7±4.3 ms at 7T (six subjects).

    pVessel.x0          = 0;   % [mm] vessel center x           | default 0
    pVessel.y0          = 0;   % [mm] vessel center y           | default 0
    pVessel.ID          = 1;   % [mm] vessel inner diameter     | default 1
    pVessel.PD          = 0;   % [mm] plug flow center diameter | default 0
    pVessel.WT          = 0;   % [mm] vessel wall thickness     | default 0
    pVessel.profile     = 'parabolic1';   %                                | default 'parabolic1'
    pVessel.vMax        = 20;   % [cm/s] maximum cross-sectional velocity | default 20
    pVessel.vMean       = 10;   % [cm/s] cross-sectional mean velocity    | default 10
    pVessel.vFlow       = [];   % [ml/min] vessel flow | default ?
    pVessel.S.lumenLami = [0.05 0.3];  % [{0,1}au] 1 values->no inflow enhancement, 2 values->linear inflow enhancement from first to second value | default blood at 7T (requires acquistion parameters)
    pVessel.S.lumenPlug = pVessel.S.lumenLami;  % same as lumenLami
    pVessel.S.wall      = 0;  %           | default 0
    pVessel.S.surround  = 0.07;  % [{0,1}au] | default blood at 7T (requires acquistion parameters)
end


    pSim

% Define simulation grid
FOVx = pSim.FOVx;
FOVy = pSim.FOVy;
nVOXx = ceil(FOVx./pSim.voxSizeX);
nVOXy = ceil(FOVy./pSim.voxSizeY);
[xGrid, yGrid,   ~, pSim.dx, pSim.dy, pSim.nSpin, pSim.spinDensity] = setGrid(FOVx, FOVy, nVOXx, nVOXy, pSim.nSpin);


% Define velocity encoding
pVenc.mVenc; % velocity encoding bandwidth [cm/s]
pVenc.dVenc; % velocity encoding resolution [cm/s]
if isfield(pVenc,'vencList') && ~isempty(pVenc.vencList)
    venc = pVenc.vencList;
else
    if pVenc.dVenc==inf
        venc = pVenc.mVenc;
    else
        M1   = linspace(0, vencToM1(pVenc.dVenc), round(vencToM1(pVenc.dVenc)/vencToM1(pVenc.mVenc))+1);
        venc = M1toVenc(M1);
    end
    pVenc.vencList = venc;
end



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
    [If(:,:,iVenc),pSim.voxIdx] = getVoxelSignal(spinMap,pSim       ,mask.lumen   );
    [Is(:,:,iVenc),          ~] = getVoxelSignal(spinMap,pSim.voxIdx,mask.surround);
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

    if isfield(pSim,'voxRndOffsetX') && isfield(pSim,'voxRndOffsetY') && (~isempty(pSim.voxRndOffsetX) || ~isempty(pSim.voxRndOffsetY))
        x0 = pSim.voxRndOffsetX;
        y0 = pSim.voxRndOffsetY;
    else
        xSpinRange = (pSim.voxSizeX./pSim.dx-1)/2 .*[-1 1];
        ySpinRange = (pSim.voxSizeY./pSim.dy-1)/2 .*[-1 1];
        x0 = randi(xSpinRange,[1 pSim.voxRndOffset]).*pSim.dx;
        y0 = randi(ySpinRange,[1 pSim.voxRndOffset]).*pSim.dy;
    end

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
            [If(:,:,iVenc,iPos),~] = getVoxelSignal(spinMap,pSim.voxIdx,mask.lumen   );
            [Is(:,:,iVenc,iPos),~] = getVoxelSignal(spinMap,pSim.voxIdx,mask.surround);
        end
    end

    res.pSim.voxRndOffsetX = x0;
    res.pSim.voxRndOffsetY = y0;
    res.rndOffset.If = If;
    res.rndOffset.Is = Is;
    if verbose; disp('Vessel at random positions. Done.'); end
end
