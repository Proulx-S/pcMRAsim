function res = runSim(pVessel, pSim, pMri, verbose)
if ~exist('verbose','var') || isempty(verbose); verbose = true; end

%% Default parameters
% Vessel simulation parameters
if ~exist('pVessel','var') || isempty(pVessel)
    % vessel geometry and flow
    pVessel.ID          = 1;   % [mm]     vessel inner diameter
    pVessel.PD          = 0;   % [mm]     plug flow center diameter
    pVessel.WT          = 0;   % [mm]     vessel wall thickness
    pVessel.profile     = 'parabolic1'; % flow profile
    pVessel.vMax        = 20;  % [cm/s]   maximum cross-sectional (through-slice) velocity
    pVessel.vMean       = 10;  % [cm/s]   mean    cross-sectional (through-slice) velocity
    pVessel.vFlow       = [];  % [ml/min] blood flow
    % mr signal intensities -- leave empty for a determination based on relaxation and acquisition parameters
    pVessel.S.lumen     = [];  % [MR signal {0,1}] from the vessel lumen    compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.wall      = 0;   % [MR signal {0,1}] from the vessel wall     compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.surround  = [];  % [MR signal {0,1}] from the static surround compartment if it filled the whole voxel | determined from pMri if unspecified
end
% Spin simulation parameters
if ~exist('pSim','var') || isempty(pSim)
    pSim.voxFE       = 1;          % [mm]     voxel  size in FE direction
    pSim.voxPE       = pSim.voxFE; % [mm]     voxel  size in PE direction
    pSim.matFE       = 3;          % [voxels] matrix size in FE direction (must be odd)
    pSim.matPE       = pSim.matFE; % [voxels] matrix size in PE direction (must be odd)
    pSim.nSpin       = (2^8)^2;    % [n] spins per voxel
    % randomization of vessel position relative to center voxel
    pSim.monteCarloN = 0;          % [n] number of bootstrap object-to-grid random shifts
end
% MR parameters
if ~exist('pMri','var') || isempty(pMri)
    % imaging
    pMri.sliceThickness     = 1;      % [mm]
    pMri.TR                 = 0.05;   % [s]   RF repetition time (alpha TR)
    pMri.TE                 = 0.008;  % [s]   echo time
    pMri.FA                 = 40;     % [deg]
    % velocity encoding
    pMri.venc.method = 'PCmono'; % 'FVEmono' | 'FVEbipo' | 'PCmono' | 'PCbipo'
    switch pMri.venc.method
        case 'FVEmono'    % monopolar fourier velocity encoding
            pMri.venc.vencRes       = 2;                % [cm/s]    velocity spectrum resolution (minimum velocity)
            pMri.venc.vencMax       = 50;               % [cm/s]    velocity spectrum span       (maximum velocity)
            m1List = linspace(0, vencToM1(pMri.venc.vencRes), round(vencToM1(pMri.venc.vencRes)/vencToM1(pMri.venc.vencMax))+1);
            pMri.venc.m1List        = m1List;           % [T*s^2/m] list of velocity encoding gradient first moments
            pMri.venc.vencList      = M1toVenc(m1List); % [cm/s]    list of velocity encoding values
        case 'FVEbipo'    % bipolar fourier velocity encoding
            pMri.venc.vencRes       = 2;                % [cm/s]    velocity spectrum resolution (minimum velocity)
            pMri.venc.vencMax       = 50;               % [cm/s]    velocity spectrum span       (maximum velocity)
            m1List = linspace(0, vencToM1(pMri.venc.vencRes), round(vencToM1(pMri.venc.vencRes)/vencToM1(pMri.venc.vencMax))+1);
            m1List = cat(2,-flip(m1List(2:end)),m1List);
            pMri.venc.m1List        = m1List;           % [T*s^2/m] list of velocity encoding gradient first moments
            pMri.venc.vencList      = M1toVenc(m1List); % [cm/s]    list of velocity encoding values
        case 'PCmono' % monopolar phase-Contrast velocity encoding
            pMri.venc.vencList      = [inf 8];          % [cm/s] list of velocity encoding values
            pMri.venc.m1List        = vencToM1(pMri.venc.vencList); % [T*s^2/m] list of velocity encoding gradient first moments
            pMri.venc.vencRes       = [];               % not used
            pMri.venc.vencMax       = [];               % not used
        case 'PCbipo' % bipolar phase-Contrast velocity encoding
            pMri.venc.vencList      = [-8 8];           % [cm/s] list of velocity encoding values
            pMri.venc.m1List        = vencToM1(pMri.venc.vencList); % [T*s^2/m] list of velocity encoding gradient first moments
            pMri.venc.vencRes       = [];               % not used
            pMri.venc.vencMax       = [];               % not used
        otherwise
            error('Invalid velocity encoding method: %s', pMri.venc.method);
    end
    % relaxation
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
end

if nargin == 0
    res.pVessel = pVessel;
    res.pSim    = pSim;
    res.pMri    = pMri;
    return;
end

% Define simulation grid
[pSim.gridFE, pSim.gridPE, pSim.gridVoxIdx, pSim.dFE, pSim.dPE, pSim.nSpin] = setGrid(pSim.voxFE, pSim.voxPE, pSim.matFE, pSim.matPE, pSim.nSpin);

% Define velocity encoding
switch pMri.venc.method
    case 'FVEmono'
        m1List = linspace(0, vencToM1(pMri.venc.vencRes), round(vencToM1(pMri.venc.vencRes)/vencToM1(pMri.venc.vencMax))+1);
        pMri.venc.m1List  = m1List;
        pMri.venc.vencList = M1toVenc(m1List);
    case 'FVEbipo'
        if ~isfield(pMri.venc,'m1List') || isempty(pMri.venc.m1List)
            pMri.venc.m1List   = linspace(0, vencToM1(pMri.venc.vencRes), round(vencToM1(pMri.venc.vencRes)/vencToM1(pMri.venc.vencMax))+1);
            pMri.venc.m1List   = cat(2,-flip(pMri.venc.m1List(2:end)),pMri.venc.m1List); % [T*s^2/m]
        end
        if ~isfield(pMri.venc,'vencList') || isempty(pMri.venc.vencList)
            pMri.venc.vencList = M1toVenc(pMri.venc.m1List); % [cm/s]    list of velocity encoding values
        end

    case 'PCmono'
    case 'PCbipo'
    otherwise
        error('Invalid velocity encoding method: %s', pMri.venc.method);
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% NEED TO ADAPT simVesselSpins.m FOR MULTIPLE VENC VALUES
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Simulate with vessel centered on center voxel
[res.magMap,res.vMap,res.pVessel,res.pSim,res.pMri] = simVesselSpins(pVessel, pSim, pMri);
res.spinMap = res.magMap.*exp(1i*vel2phase(res.vMap, res.pMri.venc.vencList(2)));
res.I  = sum(res.spinMap(res.pSim.gridVoxIdx==0                            )); % total signal
res.If = sum(res.spinMap(res.pSim.gridVoxIdx==0 & res.pVessel.mask.lumen   )); % lumen signal
res.Is = sum(res.spinMap(res.pSim.gridVoxIdx==0 & res.pVessel.mask.surround)); % surround signal

% Simulate with random position of the vessel center within the center voxel
%this will use values precomputed from above and just move the vessel around on each monte carlo iteration
if pSim.monteCarloN > 0
    dbstack; error('Not implemented, just copy-pasted from older code.');
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
