function res = runSim(pVessel, pSim, pMri, verbose, light, earlyStop)
if ~exist('verbose'  ,'var') || isempty(verbose  ); verbose   = true;  end
if ~exist('light'    ,'var') || isempty(light    ); light     = true;  end
if ~exist('earlyStop','var') || isempty(earlyStop); earlyStop = false; end

% Recovery mode: res struct passed as first arg — re-run simVesselSpins to restore magMap/vMap
if isstruct(pVessel) && isfield(pVessel, 'pVessel')
    resOld  = pVessel;
    resNew  = runSim(resOld.pVessel, resOld.pSim, resOld.pMri, verbose, light, true);
    resOld.magMap = resNew.magMap;
    resOld.vMap   = resNew.vMap;
    res = resOld;
    return;
end

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
    pVessel.vMean       = 5;  % [cm/s]   mean    cross-sectional (through-slice) velocity
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
    if exist('pMri','var') && iscell(pMri)
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
                    % pMri.relax.blood.T1     = 2.7  ; % [s]  new value = 2.900s
                    % pMri.relax.blood.T2star = 10e-3; % [s]  new value = 0.005s
                    % pMri.relax.GM.T1        = 2.3  ; % [s]  new value = 2.100s
                    % pMri.relax.GM.T2star    = 15e-3; % [s]

                    % just take it from Bates et al. 2023 (doi:10.1007/s10334-023-01081-3)
                    pMri.relax.blood.T1     = 3   ; % [s]
                    pMri.relax.blood.T2star = 5e-3; % [s]  venous T2 = 6ms, here we eyeball a shorter T2* midway between venous and arterial
                    pMri.relax.GM.T1        = 2.3  ; % [s]  
                    pMri.relax.GM.T2star    = 16e-3; % [s]



                otherwise
                    error('Invalid species: %s', pMri.species);
            end
            % Mouse at 14T: no direct measurement. Linear field dependence (Dobre et al., Magn Reson Imag 25(5):733–735, 2007):
            %   T1(ms) = 129*B0 + 1167 (1.5–9.4 T); extrapolation at 14T → 2.97 s. Reduced for higher mouse Hct (~0.48) vs human → ~2.7 s.
            %   Using 2.7 s as nominal mouse blood T1 at 14T.
            % Mouse at 14T: no direct 14T. Cortical (isocortex) at 11.7T ~2.04 s in vivo (Kumar et al., NMR Biomed 39(1):e70187, 2026, doi:10.1002/nbm.70187).
            %   T1 increases with B0; extrapolation 11.7T→14T → ~2.3 s. Using 2.3 s as nominal mouse cortical GM T1 at 14T.
            % Mouse at 14T: venous T2 (not T2*) at 11.7T 26.9±1.7 ms normoxia (Wei et al., MRM 80:521–528, 2018, doi:10.1002/mrm.27046).
            %   T2* < T2 due to susceptibility; T2* shortens with B0. Estimate venous T2* at 14T ~10 ms.
            % Mouse at 14T: no direct 14T. R2* increases ~linearly with B0; at 17.6T mouse brain T2* measured (Kara et al., MRM 70:985–993, 2013).
            %   Extrapolation 7T→14T: T2* scales roughly as 1/B0 → cortical GM at 14T ~15 ms. Using 15 ms as nominal.


            
            if 0
            % --------BLOOD T1--------
            HctA = median([0.45 0.42 0.41 0.43 0.43 0.53]); % 0.43 (from Dobre)
            HctB = 0.48; %(for mouse)
            % ----doi:10.1002/mrm.20178
            % --arterial: ~69ms difference
            R1a=0.52.*HctA+0.38;
            R1b=0.52.*HctB+0.38;
            1/R1a-1/R1b
            % --venous: ~96ms difference
            R1a=0.83.*HctA+0.28;
            R1b=0.83.*HctB+0.28;
            (1/R1a-1/R1b)
            % ----doi:10.1002/mrm.24547
            % --arterial: ~85ms difference
            R1a = (0.305+0.201*HctA) / (0.95-0.25*HctA);
            R1b = (0.305+0.201*HctB) / (0.95-0.25*HctB);
            (1/R1a-1/R1b)
            % --venous: ~91ms difference
            R1a = (0.305+0.275*HctA) / (0.95-0.25*HctA);
            R1b = (0.305+0.275*HctB) / (0.95-0.25*HctB);
            (1/R1a-1/R1b)
            end

            
            
            if 0
            % --------BLOOD T2--------
            B0=[1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	3	3	3	3	3	3	3	3	3	7	7	1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	1.5	2.35	3	3	3	3	3	3	4.7	4.7	4.7	4.7	7	7	7	7]; %T
            T2=[172	131	209	114	172	242	190	253	61.9	68.4	155	169	72.4	67.7	80.1	60.9	80.3	19.5	27.4	212	201	159	142	131	98.1	88.6	86.9	156	145	174	245	98.9	84	72.4	66.7	63.2	72.4	63.2	12.2	110	51.5	9.1	19.4	27	24.1	22]; %ms
            [B0,b] = sort(B0);
            T2 = T2(b);
            censorIdx = 1./T2>0.06;
            figure
            plot(B0(~censorIdx),T2(~censorIdx),'o','MarkerEdgeColor','k','MarkerFaceColor','w'); hold on;
            plot(B0(censorIdx),T2(censorIdx),'o','MarkerEdgeColor','k','MarkerFaceColor','r');
            fT2b     = fit(B0(:),1./T2(:),fittype('p1*x+p2*x^2','coefficients',{'p1','p2'}));
            fT2b_cen = fit(B0(~censorIdx)',1./T2(~censorIdx)',fittype('p1*x+p2*x^2','coefficients',{'p1','p2'}));
            B0fit = linspace(0,14,200);
            plot(B0fit,1./fT2b(B0fit),'r-');
            plot(B0fit,1./fT2b_cen(B0fit),'w-');
            plot(14,1./fT2b(14),'rx','MarkerSize',10,'LineWidth',2);
            plot(14,1./fT2b_cen(14),'wx','MarkerSize',10,'LineWidth',2);
            xline(14,'w--');
            plot(11.7,32.3,'cx','MarkerSize',10,'LineWidth',2);
            grid on
            legend('data','censored data','fit','fit after censoring',sprintf('fit @14T: T2=%.1f ms',1/fT2b(14)),sprintf('censored @14T: T2=%.1f ms',1/fT2b_cen(14)),'B0 = 14T',sprintf('Wei et al. 2018: T2=%.1f ms',32.3),'Location','northeast');
            end



            if 0
            % --------GM T1--------
            B0=[0.5	1.5	1.5	1.5	1.5	2	2	3	3	4	4	4	7	7	9.4	9.4	1.5	2	4	4	4	4	4	4.7	7	9.4	9.4	9.4	9.4	9.4	11.7	11.7	11.7	11.7	17.2	17.6]; %T
            T1=[650	1050	1100	1175	1200	1200	1250	1500	1550	1250	1300	1325	1875	1900	1875	1950	1000	1200	875	950	1050	1250	1300	1600	1600	1575	1800	1850	2125	2150	1575	1600	1850	2050	2600	2000]; %ms
            kumarB0 = 11.7;
            kumarT1 = 2036;
            [B0,b] = sort(B0);
            T1 = T1(b);
            fVen = fit(B0(:),T1(:),fittype('p1*x^p2','coefficients',{'p1','p2'}));
            fRooney = fVen; fRooney.p1 = 0.857*1000; fRooney.p2 = 0.376;
            fVenAdjKumar = fVen; fVenAdjKumar.p1 = fVen.p1 * kumarT1/fVen(kumarB0);
            B0fit  = linspace(0,max(B0),200);
            B0fit2 = linspace(kumarB0,14,10);
            figure;
            plot(B0,T1,'o','MarkerEdgeColor','k','MarkerFaceColor','w'); hold on;
            plot(B0fit,fRooney(B0fit),'r-');
            plot(B0fit,fVen(B0fit),'w-');
            % plot(kumarB0,kumarT1,'cx','MarkerSize',10,'LineWidth',2);
            plot(B0fit2,fVenAdjKumar(B0fit2),'c-');
            T1_Rooney14 = fRooney(14); T1_Ven14 = fVen(14); T1_KumarAdj14 = fVenAdjKumar(14);
            xline(14,'w--');
            yline(T1_Rooney14,'r--');
            yline(T1_Ven14,'w--');
            plot(kumarB0,kumarT1,'cx','MarkerSize',10,'LineWidth',2);
            plot(14,T1_KumarAdj14,'cx','MarkerSize',10,'LineWidth',2);
            grid on; axis square;
            ylabel('T1 [ms]'); legend('data','Rooney  et al. 2007 model','van de Ven et al. 2007 fit','Kumar adjusted van de Ven fit',sprintf('B0 = 14T'),sprintf('Rooney T1: %.0f ms',T1_Rooney14),sprintf('van de Ven T1: %.0f ms',T1_Ven14),sprintf('Kumar et al. 2026: %.0f ms @%.1fT',kumarT1,kumarB0),sprintf('Kumar adj. @14T: %.0f ms',T1_KumarAdj14),'Location','southeast');
            end

            
            


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
    case 'PCmono'
        pMri.venc.vencList = pMri.venc.vencList(:);
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
        pMri.venc.m1List = cat(2,pMri.venc.m1List,zeros(size(pMri.venc.m1List))); % second line for references (M1=0 in the monopolar case)
    case 'PCbipo'
        pMri.venc.vencList = pMri.venc.vencList(:);
        pMri.venc.m1List = vencToM1(pMri.venc.vencList); % [T*s^2/m]
        pMri.venc.m1List = cat(2,pMri.venc.m1List,-pMri.venc.m1List)./2; % second line for references (-M1 in the bipolar case) and divide by 2 for bipolar encoding
    otherwise
        error('Invalid velocity encoding method: %s', pMri.venc.method);
end

% Simulate with vessel centered on center voxel
[res.magMap,res.vMap,res.pVessel,res.pSim,res.pMri] = simVesselSpins(pVessel, pSim, pMri);
if earlyStop; return; end

% Simulate spin map
m1 = permute(res.pMri.venc.m1List,[3 4 5 6 1 2 7 8 9 10 11 12 13 14 15 16]);
res.spinMap = res.magMap.*exp(1i*vel2phase(res.vMap, m1));
spinMap = permute(res.spinMap,[5 6 7 8 9 10 11 12 13 14 15 16 1 2 3 4]);



%% Center voxel averaging
I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0                            ),13); % total signal
If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0 & res.pVessel.mask.lumen   ),13); % lumen signal
Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,res.pSim.gridVoxIdx==0 & res.pVessel.mask.surround),13); % surround signal
res.I  = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
res.If = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
res.Is = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);

dimList = {'FE' 'PE' 'SL' 't' 'M1' 'M1ref'};
res.info = strjoin(dimList,' x ');

% Simulate with random position of the vessel center within the center voxel
%this will use values precomputed from above and just move the vessel around on each monte carlo iteration
if pSim.monteCarloN > 0

    res.I  = cat(7,res.I ,nan([size(res.I ,1:6) pSim.monteCarloN]));
    res.If = cat(7,res.If,nan([size(res.If,1:6) pSim.monteCarloN]));
    res.Is = cat(7,res.Is,nan([size(res.Is,1:6) pSim.monteCarloN]));
    res.info = strjoin({res.info 'mntCrls'},' x ');

    for iMntCrl = 1:pSim.monteCarloN

        % shift the grid voxels by the monte carlo shift, relative to the spin map
        gridVoxIdx = circshift(res.pSim.gridVoxIdx, [res.pSim.monteCarloShiftFE(iMntCrl) res.pSim.monteCarloShiftPE(iMntCrl)] );

        I  = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0                            ),13); % total signal
        If = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0 & res.pVessel.mask.lumen   ),13); % lumen signal
        Is = sum(spinMap(:,:,:,:,:,:,:,:,:,:,:,:,gridVoxIdx==0 & res.pVessel.mask.surround),13); % surround signal

        res.I( :,:,:,:,:,:,iMntCrl+1) = permute(I,  [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.If(:,:,:,:,:,:,iMntCrl+1) = permute(If, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
        res.Is(:,:,:,:,:,:,iMntCrl+1) = permute(Is, [13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12]);
    end
    if verbose; disp('Vessel at random positions. Done.'); end
end



%% Subtract ref phase
switch pMri.venc.method
    case {'FVEmono','FVEbipo'}
    case {'PCmono' 'PCbipo'}
        res.I  = res.I  ./ exp(1i*angle(res.I( :,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
        res.If = res.If ./ exp(1i*angle(res.If(:,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
        res.Is = res.Is ./ exp(1i*angle(res.Is(:,:,:,:,:,end,:,:,:,:,:,:,:,:,:,:)));
        res.info2 = 'ref phase subtracted';
end



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
res.pSim.gridFE     = single(res.pSim.gridFE   );
res.pSim.gridPE     = single(res.pSim.gridPE   );
res.pSim.gridVoxIdx = uint8(res.pSim.gridVoxIdx);
