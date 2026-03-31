function res = runSim(pVessel, pSim, pMri, verbose)
if ~exist('verbose','var') || isempty(verbose); verbose = true; end

    % when pVessel.S.lumen and pVessel.S.surround are empty, they are determined from the relaxation and acquisition parameters
    % when pVessel.profile is numeric, it is used as the velocity profile -- this allows to specify and arbitrary velocity distribution within the whole ROI (not the center voxel). Must be of length pSim.nSpin.

%% Default parameters
% Vessel simulation parameters
if ~exist('pVessel','var') || isempty(pVessel)
    % vessel geometry and flow
    pVessel.ID          = 1;   % [mm]     vessel inner diameter
    pVessel.PD          = 0;   % [mm]     plug flow center diameter
    pVessel.WT          = 0;   % [mm]     vessel wall thickness
    pVessel.profile     = 'parabolic1'; % flow profile: 'parabolic' | 'parabolic1' | 'plug' | 'plugFlow'
    pVessel.vMean       = 10;  % [cm/s]   mean    cross-sectional (through-slice) velocity
    pVessel.vMax        = [];  % [cm/s]   maximum cross-sectional (through-slice) velocity
    pVessel.vFlow       = [];  % [ml/min] blood flow
    % mr signal intensities -- leave empty for a determination based on relaxation and acquisition parameters
    pVessel.S.lumen     = [];  % [MR signal {0,1}] from the vessel lumen    compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.wall      = 0;   % [MR signal {0,1}] from the vessel wall     compartment if it filled the whole voxel | determined from pMri if unspecified
    pVessel.S.surround  = [];  % [MR signal {0,1}] from the static surround compartment if it filled the whole voxel | determined from pMri if unspecified
end
% Spin simulation parameters
if ~exist('pSim','var') || isempty(pSim)
    pSim.fovFE       = 3;          % [mm]     field of view in FE direction
    pSim.fovPE       = pSim.fovFE; % [mm]     field of view in PE direction
    pSim.matFE       = 3;          % [voxels] matrix size in FE direction (must be odd)
    pSim.matPE       = pSim.matFE; % [voxels] matrix size in PE direction (must be odd)
    pSim.nSpin       = (2^8+1)^2;    % [n] spins per voxel
    % randomization of vessel position relative to center voxel
    pSim.monteCarloN = 0;          % [n] number of bootstrap object-to-grid random shifts
end
% MR parameters
if ~exist('pMri','var') || iscell(pMri) || isempty(pMri)
    % imaging
    if iscell(pMri)
        vencMethodList = {'FVEmono','FVEbipo','PCmono','PCbipo'};
        fieldStrengthList = {'7t','14t'};
        vencMethod    = pMri{ismember(pMri, vencMethodList)};
        fieldStrength = pMri{ismember(pMri, fieldStrengthList)};
        clear pMri;
    else
        fieldStrength = '7t';
        vencMethod    = 'FVEbipo';
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
        case 'FVEmono'    % monopolar fourier velocity encoding
            pMri.venc.FVEres       = 2;                % [cm/s]    velocity spectrum resolution (minimum velocity)
            pMri.venc.FVEbw        = 50;               % [cm/s]    velocity spectrum span       (maximum velocity)
            [pMri.venc.vencList, pMri.venc.m1List, pMri.venc.FVEvel, pMri.venc.Ns, pMri.venc.vencMin, pMri.venc.vencMax] = getFVE(pMri.venc.FVEres, pMri.venc.FVEbw, pMri.venc.method);
            pMri.venc.vencList;                         % [cm/s]    list of velocity encoding values
            pMri.venc.FVEvel;                           % [cm/s]    velocity spectrum "frequency" axis
            pMri.venc.vencMin;                          % [cm/s]    minimum venc value
            pMri.venc.vencMax;                          % [cm/s]    maximum venc value
        case 'FVEbipo'    % bipolar fourier velocity encoding
            pMri.venc.FVEres       = 2;                % [cm/s]    velocity spectrum resolution (minimum velocity)
            pMri.venc.FVEbw        = 50;               % [cm/s]    velocity spectrum span       (maximum velocity)
            [pMri.venc.vencList, pMri.venc.m1List, pMri.venc.FVEvel, pMri.venc.Ns, pMri.venc.vencMin, pMri.venc.vencMax] = getFVE(pMri.venc.FVEres, pMri.venc.FVEbw, pMri.venc.method);
            pMri.venc.vencList;                         % [cm/s]    list of velocity encoding values
            pMri.venc.FVEvel;                           % [cm/s]    velocity spectrum "frequency" axis
            pMri.venc.vencMin;                          % [cm/s]    minimum venc value
            pMri.venc.vencMax;                          % [cm/s]    maximum venc value
        case 'PCmono' % monopolar phase-Contrast velocity encoding
            pMri.venc.FVEres = 0; % not used
            pMri.venc.FVEbw  = 0; % not used
            pMri.venc.vencList = [inf 8]';         % [cm/s] list of velocity encoding values
            pMri.venc.FVEres  = []; % not used
            pMri.venc.FVEbw   = []; % not used
            pMri.venc.FVEvel  = []; % not used
            pMri.venc.vencMin = []; % not used
            pMri.venc.vencMax = []; % not used
        case 'PCbipo' % bipolar phase-Contrast velocity encoding
            pMri.venc.FVEres  = 0; % not used
            pMri.venc.FVEbw   = 0; % not used
            pMri.venc.vencList = [-8 8]';           % [cm/s] list of velocity encoding values
            pMri.venc.FVEres  = []; % not used
            pMri.venc.FVEbw   = []; % not used
            pMri.venc.FVEvel  = []; % not used
            pMri.venc.vencMin = []; % not used
            pMri.venc.vencMax = []; % not used
        otherwise
            error('Invalid velocity encoding method: %s', pMri.venc.method);
    end
    % relaxation
    switch pMri.fieldStrength
        case 7
            switch pMri.species
                case 'human'
                    pMri.relax.blood.T1     = 2.58   ;   % [s]
                    pMri.relax.blood.T2star = 10e-3  ;   % [s]
                    pMri.relax.GM.T1        = 1.939  ;  % [s]
                    pMri.relax.GM.T2star    = 32.9e-3; % [s]
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
        otherwise
            error('Invalid field strength: %s', pMri.fieldStrength);
    end

    
end

if nargin == 0
    res.pVessel = pVessel;
    res.pSim    = pSim;
    res.pMri    = pMri;
    return;
end

% Define simulation grid
[pSim.gridFE, pSim.gridPE, pSim.gridVoxIdx, pSim.dFE, pSim.dPE, pSim.nSpin] = setGrid(pSim.fovFE, pSim.fovPE, pSim.matFE, pSim.matPE, pSim.nSpin);

% Define velocity encoding
switch pMri.venc.method
    case {'FVEmono','FVEbipo'}
        [pMri.venc.vencList, pMri.venc.m1List, pMri.venc.FVEvel, pMri.venc.Ns, pMri.venc.vencMin, pMri.venc.vencMax] = getFVE(pMri.venc.FVEres, pMri.venc.FVEbw, pMri.venc.method);
        pMri.venc.m1List; % [T*s^2/m]
    case {'PCmono' 'PCbipo'}
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
    otherwise
        error('Invalid velocity encoding method: %s', pMri.venc.method);
end

% Simulate with vessel centered on center voxel
[res.magMap,res.vMap,res.pVessel,res.pSim,res.pMri] = simVesselSpins(pVessel, pSim, pMri);

% Simulate spins and center voxel averaging
res.spinMap = res.magMap.*exp(1i*vel2phase(res.vMap, res.pMri.venc.vencList));
spinMap = permute(res.spinMap,[5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4]);
I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0                            ),13); % total signal
If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0 & res.pVessel.mask.lumen   ),13); % lumen signal
Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0 & res.pVessel.mask.surround),13); % surround signal
res.I  = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
res.If = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
res.Is = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);

dimList = {'FE' 'PE' 'SL' 't' 'venc'};
res.info = strjoin(dimList,' x ');



% Simulate with random position of the vessel center within the center voxel
%this will use values precomputed from above and just move the vessel around on each monte carlo iteration
if pSim.monteCarloN > 0

    res.I  = cat(6,res.I ,nan([size(res.I ,1:5) pSim.monteCarloN]));
    res.If = cat(6,res.If,nan([size(res.If,1:5) pSim.monteCarloN]));
    res.Is = cat(6,res.Is,nan([size(res.Is,1:5) pSim.monteCarloN]));
    res.info = strjoin({res.info 'mntCrls'},' x ');


    nSpinFE = max(sum(res.pSim.gridVoxIdx==0,2));
    shiftFE = (1:nSpinFE)-nSpinFE/2-0.5;
    nSpinPE = max(sum(res.pSim.gridVoxIdx==0,1));
    shiftPE = (1:nSpinPE)-nSpinPE/2-0.5;    

    for iMntCrl = 1:pSim.monteCarloN
        gridVoxIdx = circshift(res.pSim.gridVoxIdx, [shiftFE(randperm(length(shiftFE),1)) shiftPE(randperm(length(shiftPE),1))] );

        I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0                            ),13); % total signal
        If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0 & res.pVessel.mask.lumen   ),13); % lumen signal
        Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0 & res.pVessel.mask.surround),13); % surround signal

        res.I( :,:,:,:,:,iMntCrl+1) = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.If(:,:,:,:,:,iMntCrl+1) = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.Is(:,:,:,:,:,iMntCrl+1) = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
    end
    if verbose; disp('Vessel at random positions. Done.'); end
end
