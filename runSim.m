function res = runSim(pVessel, pSim, pMri, verbose, light, earlyStop)
if ~exist('verbose'  ,'var') || isempty(verbose  ); verbose   = true;  end
if ~exist('light'    ,'var') || isempty(light    ); light     = true;  end
if ~exist('earlyStop','var') || isempty(earlyStop); earlyStop = false; end

% Recovery mode: res struct passed as first arg — re-run simVesselSpins to restore magMap/vMap
if nargin~=0 && isstruct(pVessel) && isfield(pVessel, 'pVessel')
    resOld  = pVessel;
    resNew  = runSim(resOld.pVessel, resOld.pSim, resOld.pMri, verbose, light, true);
    resOld.magMap = resNew.magMap;
    resOld.vMap   = resNew.vMap;
    res = resOld;
    return;
end

    % when pVessel.S.lumen and pVessel.S.surround are empty, they are determined from the relaxation and acquisition parameters
    % when pVessel.profile is numeric, it is used as the velocity profile -- this allows to specify and arbitrary velocity distribution within the whole ROI (not the center voxel). Must be of length pSim.spinGrid.matFE*pSim.spinGrid.matPE.

%% Default parameters
% Vessel simulation parameters
if ~exist('pVessel','var') || isempty(pVessel)
    % vessel geometry and flow
    pVessel.ID          = 5;   % [mm]     vessel inner diameter
    pVessel.PD          = 0;   % [mm]     plug flow center diameter
    pVessel.WT          = 0;   % [mm]     vessel wall thickness
    pVessel.profile     = 'parabolic1'; % flow profile: 'parabolic' | 'parabolic1' | 'plug' | 'plugFlow'
    pVessel.vMean       = 1.5;  % [cm/s]   mean    cross-sectional (through-slice) velocity
    pVessel.vMax        = [];  % [cm/s]   maximum cross-sectional (through-slice) velocity
    pVessel.vFlow       = [];  % [ml/min] blood flow
    % mr signal intensities -- leave empty for a determination based on relaxation and acquisition parameters
    pVessel.S.lumen     = [];  % [MR signal {0,1}] from the vessel lumen    compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.wall      = 0;   % [MR signal {0,1}] from the vessel wall     compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.surround  = [];  % [MR signal {0,1}] from the static surround compartment if it filled the whole voxel | determined from pMri if unspecified
    % pVessel.posFE       = 0;   % [mm]     vessel center FE offset from grid center
    % pVessel.posPE       = 0;   % [mm]     vessel center PE offset from grid center
    pVessel.pos   = [0 0 0]; % [mm]     vessel center offset from grid center
    pVessel.n_hat = [0 0 1]; % [a.u.] vessel center offset from grid center
end
% Spin simulation parameters
if ~exist('pSim','var') || isempty(pSim)
    pSim.voxGrid.fovFE    = 3;          % [mm]     field of view in FE  direction
    pSim.voxGrid.fovPE    = 3;          % [mm]     field of view in PE  direction
    pSim.voxGrid.fovSLC   = 2;          % [mm]     field of view in SLC direction
    pSim.voxGrid.matFE    = 3;          % [voxels] matrix size in FE  direction
    pSim.voxGrid.matPE    = 3;          % [voxels] matrix size in PE  direction
    pSim.voxGrid.matSLC   = 1;          % [voxels] matrix size in SLC direction
    pSim.isochromatPerVoxPerDirMin = 2; % [n]      minimum number of isochromat in each direction in a voxel
    pSim.isochromatPerVox = [];   % [n]      target total spin count across grid
    pSim.gridMode      = 'pseudoVoxel'; % 'pseudoVoxel' | 'centerVox' (centerVox requires odd matFE/matPE)
    pSim.monteCarloN   = 0;          % [n]      Monte Carlo vessel-position shifts (centerVox only)
    pSim.mriScale      = 1;      %
end
% MR parameters
if ~exist('pMri','var') || iscell(pMri) || isempty(pMri)
    % imaging
    if exist('pMri','var') && iscell(pMri)
        vencMethodList = {'FVEmono','FVEbipo','PCmono','PCbipo','velEnc'};
        fieldStrengthList = {'7t','14t'};
        vencMethod    = pMri{ismember(pMri, vencMethodList)};
        fieldStrength = pMri{ismember(pMri, fieldStrengthList)};
        clear pMri;
    else
        fieldStrength = '7t';
        vencMethod    = 'FVEmono';
    end
    switch lower(fieldStrength)
        case '7t'
            pMri.fieldStrength = 7;
            pMri.species       = 'human';
        case '14t'
            pMri.fieldStrength = 14;
            pMri.species       = 'mouse';
        otherwise
            error('Invalid field strength: %s', fieldStrength);
    end
    clear fieldStrength;
    pMri.sliceThickness     = 1;      % [mm]
    pMri.TR                 = 0.05;   % [s]   RF repetition time (alpha TR)
    pMri.TE                 = 0.008;  % [s]   echo time
    pMri.FA                 = 40;     % [deg]
    % velocity encoding
    pMri.venc.method = vencMethod; % 'FVEmono' | 'FVEbipo' | 'PCmono' | 'PCbipo'
    clear vencMethod;
    switch pMri.venc.method
        case {'FVEmono' 'FVEbipo'}    % monopolar/bipolar fourier velocity encoding
            pMri.venc.FVEres       = 2;                % [cm/s]    velocity spectrum resolution (minimum velocity)
            pMri.venc.FVEbw        = 50;               % [cm/s]    velocity spectrum span       (maximum velocity)
            [pMri.venc.vencList, pMri.venc.m1List, pMri.venc.FVEvel, pMri.venc.Ns, pMri.venc.vencMin, pMri.venc.vencMax] = getFVE(pMri.venc.FVEres, pMri.venc.FVEbw, pMri.venc.method);
            pMri.venc.vencList;                         % [cm/s]    list of velocity encoding values
            pMri.venc.FVEvel;                           % [cm/s]    velocity spectrum "frequency" axis
            pMri.venc.vencMin;                          % [cm/s]    minimum venc value
            pMri.venc.vencMax;                          % [cm/s]    maximum venc value
        case {'PCmono' 'PCbipo'} % monopolar/bipolar (one-sided/two-sided) phase-contrast velocity encoding
            pMri.venc.FVEres   = 0;  % not used
            pMri.venc.FVEbw    = 0;  % not used
            pMri.venc.vencList = [20 40]; % [cm/s] list of velocity encoding values
            pMri.venc.FVEres   = []; % not used
            pMri.venc.FVEbw    = []; % not used
            pMri.venc.FVEvel   = []; % not used
            pMri.venc.vencMin  = []; % not used
            pMri.venc.vencMax  = []; % not used
        otherwise
            error('Invalid velocity encoding method: %s', pMri.venc.method);
    end


end
% Populate relaxation parameters from fieldStrength and species.
% Runs unconditionally so that changing fieldStrength/species on a pre-built pMri struct
% is picked up on the next runSim call without having to clear pMri.relax manually.
if isfield(pMri,'fieldStrength') && isfield(pMri,'species')
    switch pMri.fieldStrength
        case 7
            switch pMri.species
                case 'human'
                    pMri.relax.blood.T1     = 2.58   ;   % [s]
                    pMri.relax.blood.T2star = 10e-3  ;   % [s]
                    pMri.relax.GM.T1        = 1.939  ;   % [s]
                    pMri.relax.GM.T2star    = 32.9e-3;   % [s]
                otherwise
                    error('Invalid species: %s', pMri.species);
            end
            % Blood T1. Human at 7T (Rane & Gore, Magn Reson Imaging 31(3):477–479, 2013, doi:10.1016/j.mri.2012.08.008):
            %   arterial 2.29±0.10 s, venous 2.07±0.12 s in vitro (37°C); venous sagittal sinus in vivo 2.45±0.11 s.
            %   arterial in vivo 2.45/2.07*2.29 = 2.71 s
            %   mid arterio-venous in vivo (2.71+2.45)/2 = 2.58 s
            % Blood T2*. Human at 7T: venous blood ~7.4 ms (SWI/venography at 7T); arterial longer (higher oxygenation).
            %   Blood T2* is strongly oxygenation-dependent; R2* increases with deoxyhemoglobin (e.g. Qin & van Zijl, MRM 24868, 2009).
            % Gray matter cortical T1. Human at 7T (Waddell et al., MAGMA 21:121–130, 2008, doi:10.1007/s10334-008-0104-8):
            %   cortical gray matter 1939±149 ms (~1.94 s), white matter 1126±97 ms (MPRAGE, 4 subjects).
            % Gray matter cortical T2*. Human at 7T (Peters et al., Magn Reson Imaging 25:748–753, 2007, doi:10.1016/j.mri.2007.02.014):
            %   cortical gray matter 32.9±2.3 ms, white matter 27.7±4.3 ms at 7T (six subjects).
        case 14
            switch pMri.species
                case 'mouse'
                    pMri.relax.blood.T1     = 2.7  ; % [s]
                    pMri.relax.blood.T2star = 10e-3; % [s]
                    pMri.relax.GM.T1        = 2.3  ; % [s]
                    pMri.relax.GM.T2star    = 15e-3; % [s]
                otherwise
                    error('Invalid species: %s', pMri.species);
            end
            % Mouse at 14T: no direct measurement. Linear field dependence (Dobre et al., Magn Reson Imag 25(5):733–735, 2007):
            %   T1(ms) = 129*B0 + 1167 (1.5–9.4 T); extrapolation at 14T → 2.97 s. Reduced for higher mouse Hct (~0.48) vs human → ~2.7 s.
            %   Using 2.7 s as nominal mouse blood T1 at 14T.
            % Mouse at 14T: no direct 14T. Cortical (isocortex) at 11.7T ~2.04 s reported in vivo (e.g. wildtype mouse at 11.7T).
            %   T1 increases with B0; extrapolation 11.7T→14T → ~2.3 s. Using 2.3 s as nominal mouse cortical GM T1 at 14T.
            % Mouse at 14T: venous T2 (not T2*) at 11.7T 26.9±1.7 ms normoxia (Wei et al., MRM 80:521–528, 2018, doi:10.1002/mrm.27046).
            %   T2* < T2 due to susceptibility; T2* shortens with B0. Estimate venous T2* at 14T ~10 ms.
            % Mouse at 14T: no direct 14T. R2* increases ~linearly with B0; at 17.6T mouse brain T2* measured (Kara et al., MRM 70:985–993, 2013).
            %   Extrapolation 7T→14T: T2* scales roughly as 1/B0 → cortical GM at 14T ~15 ms. Using 15 ms as nominal.
        case 3
            switch pMri.species
                case 'phantom'
                    pMri.relax.blood.T1     = 3.25; % [s] TODO: use truer values
                    pMri.relax.blood.T2star = 0.25; % [s] TODO: use truer values
                    pMri.relax.GM.T1        = 1.10; % [s] TODO: use truer values
                    pMri.relax.GM.T2star    = 0.05; % [s] TODO: use truer values
                case 'human'
                    pMri.relax.blood.T1     = 1.66  ; % [s]   Dobre et al., MRM 2007
                    pMri.relax.blood.T2star = 50e-3 ; % [s]   arterial, oxygenation-dependent
                    pMri.relax.GM.T1        = 1.30  ; % [s]   3T cortical GM
                    pMri.relax.GM.T2star    = 30e-3 ; % [s]   3T cortical GM
                otherwise
                    error('Invalid species: %s', pMri.species);
            end
        otherwise
            error('Invalid field strength: %d', pMri.fieldStrength);
    end
end

% Define simulation grid
if ~isfield(pSim,'gridMode'); pSim.gridMode = 'pseudoVoxel'; end
if ~isfield(pVessel,'pos'); pVessel.pos = [0 0 0];             end
% if ~isfield(pVessel,'posPE'); pVessel.posPE = 0;             end
fov = [pSim.voxGrid.fovFE, pSim.voxGrid.fovPE, pSim.voxGrid.fovSLC];
mat = [pSim.voxGrid.matFE, pSim.voxGrid.matPE, pSim.voxGrid.matSLC];
[pSim.voxGrid, pSim.spinGrid, pSim.isochromatPerVox] = setGrid(fov, mat, pSim.isochromatPerVoxPerDirMin, pSim.gridMode);

if nargin == 0
    res.pVessel = pVessel;
    res.pSim    = pSim;
    res.pMri    = pMri;
    return;
end

% Define velocity encoding
switch pMri.venc.method
    case {'FVEmono','FVEbipo'}
        [pMri.venc.vencList, pMri.venc.m1List, pMri.venc.FVEvel, pMri.venc.Ns, pMri.venc.vencMin, pMri.venc.vencMax] = getFVE(pMri.venc.FVEres, pMri.venc.FVEbw, pMri.venc.method);
        pMri.venc.m1List; % [T*s^2/m]
    case 'PCmono'
        pMri.venc.vencList = pMri.venc.vencList(:);
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
        pMri.venc.m1List = cat(2,pMri.venc.m1List,zeros(size(pMri.venc.m1List))); % second line for references (M1=0 in the monopolar case)
    case 'PCbipo'
        pMri.venc.vencList = pMri.venc.vencList(:);
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
        pMri.venc.m1List = cat(2,pMri.venc.m1List,-pMri.venc.m1List)./2; % second line for references (-M1 in the bipolar case) and divide by 2 for bipolar encoding
    case 'velEnc'
        pMri.venc.vencList = pMri.venc.vencList(:);
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
    otherwise
        error('Invalid velocity encoding method: %s', pMri.venc.method);
end

% Simulate with vessel at specified position (default: center of grid)
pVessel.n_hat = pVessel.n_hat./norm(pVessel.n_hat);
[res.magMap,res.vMap,res.pVessel,res.pSim,res.pMri] = simVesselSpins(pVessel, pSim, pMri);
if earlyStop; return; end

% tmp = res.vMap;
% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorPE,squeeze(tmp(round(end/2),:,:))); axis image; colorbar
% imagesc(pSim.spinGrid.coorSLC,pSim.spinGrid.coorFE,squeeze(tmp(:,round(end/2),:))); axis image; colorbar
% imagesc(pSim.spinGrid.coorPE,pSim.spinGrid.coorFE ,squeeze(tmp(:,:,round(end/2)))); axis image; colorbar



% Simulate spin map
m1 = permute(res.pMri.venc.m1List,[3 4 5 6 1 2 7 8 9 10 11 12 13 14 15 16]);
if strcmp(res.pMri.venc.method,'PCbino'); dbstack; error('code that'); end
res.spinMap = res.magMap.*exp(1i*vel2phase(res.vMap, m1(:,:,:,:,:,1)));
spinMap = permute(res.spinMap,[5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4]);



%% Signal averaging
% allVox:      sum within each voxel
% centerVox:   sum within the center-voxel only (one voxel)
% pseudoVoxel: sum over all spins
switch pSim.gridMode
    case 'allVox'
        voxIdx = res.pSim.spinGrid.voxIdx;
        voxIdxList = unique(voxIdx(:));
        I          = nan(length(pMri.venc.vencList),pSim.voxGrid.matFE,pSim.voxGrid.matPE,pSim.voxGrid.matSLC);
    case 'centerVox'
        dbstack; error('code that');
        spinSel  = getVoxIdx(res.pSim.voxGrid, res.pSim.spinGrid) == 0;
    case 'pseudoVoxel'
        dbstack; error('code that');
        spinSel  = true(res.pSim.spinGrid.matFE, res.pSim.spinGrid.matPE, res.pSim.spinGrid.matSLC);
end
for i = 1:length(voxIdxList)
    I(:,i) = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:, voxIdx==voxIdxList(i)), 13);
end
res.I  = permute(I,  [2 3 4 5 1 6]);

% I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:, spinSel                               ), 13);
% If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:, spinSel & res.pVessel.mask.lumen      ), 13);
% Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:, spinSel & res.pVessel.mask.surround   ), 13);
% res.I  = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
% res.If = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
% res.Is = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);



dimList = {'FE' 'PE' 'SL' 't' 'M1' 'M1ref'};
res.info = strjoin(dimList,' x ');

% Monte Carlo: random vessel positions within center voxel (centerVox mode only)
if strcmp(pSim.gridMode,'centerVox') && pSim.monteCarloN > 0
    dbstack; error('double-check that')

    res.I  = cat(7,res.I ,nan([size(res.I ,1:6) pSim.monteCarloN]));
    res.If = cat(7,res.If,nan([size(res.If,1:6) pSim.monteCarloN]));
    res.Is = cat(7,res.Is,nan([size(res.Is,1:6) pSim.monteCarloN]));
    res.info = strjoin({res.info 'mntCrls'},' x ');

    voxIdx = getVoxIdx(res.pSim.voxGrid, res.pSim.spinGrid);
    for iMntCrl = 1:pSim.monteCarloN

        % shift the voxel map by the monte carlo shift, relative to the spin map
        voxIdxShifted = circshift(voxIdx, [res.pSim.monteCarloShiftFE(iMntCrl) res.pSim.monteCarloShiftPE(iMntCrl)]);

        I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,voxIdxShifted==0                              ),13); % total signal
        If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,voxIdxShifted==0 & res.pVessel.mask.lumen     ),13); % lumen signal
        Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,voxIdxShifted==0 & res.pVessel.mask.surround  ),13); % surround signal

        res.I( :,:,:,:,:,:,iMntCrl+1) = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.If(:,:,:,:,:,:,iMntCrl+1) = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.Is(:,:,:,:,:,:,iMntCrl+1) = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
    end
    if verbose; disp('Vessel at random positions. Done.'); end
end



% %% Subtract ref phase
% switch pMri.venc.method
%     case {'FVEmono','FVEbipo'}
%     case {'PCmono' 'PCbipo'}
%         res.I  = res.I  ./ exp(1i*angle(res.I( :,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
%         res.If = res.If ./ exp(1i*angle(res.If(:,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
%         res.Is = res.Is ./ exp(1i*angle(res.Is(:,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
%         res.info2 = 'ref phase subtracted';
% end



%% Reduce size of stored data
if light
    res.magMap  = [];
    res.vMap    = [];
    res.spinMap = [];
else
    res.magMap  = single(res.magMap );
    res.vMap    = single(res.vMap   );
    res.spinMap = single(res.spinMap);
end
res.pVessel.S.lumen    = single(res.pVessel.S.lumen   );
res.pVessel.S.wall     = single(res.pVessel.S.wall    );
res.pVessel.S.surround = single(res.pVessel.S.surround);
res.pSim.spinGrid.coorFE  = single(res.pSim.spinGrid.coorFE);
res.pSim.spinGrid.coorPE  = single(res.pSim.spinGrid.coorPE);
res.pSim.spinGrid.coorSLC = single(res.pSim.spinGrid.coorSLC);
res.pSim.voxGrid.coorFE  = single(res.pSim.voxGrid.coorFE);
res.pSim.voxGrid.coorPE  = single(res.pSim.voxGrid.coorPE);
res.pSim.voxGrid.coorSLC = single(res.pSim.voxGrid.coorSLC);
