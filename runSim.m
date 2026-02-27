function res = runSim(pVessel, pVenc, pSim, pReal, xStartRpl, yStartRpl, noAverageFlag)
    % simulate a single voxel by an evenly spaced grid of spins

    if nargin < 7
        noAverageFlag = false;
    end


    % Define vessel parameters
    % x0         = pVessel.x0;        % vessel center x [mm]
    % y0         = pVessel.y0;        % vessel center y [mm]
    ID         = pVessel.ID;        % vessel inner diameter [mm]
    WT         = pVessel.WT;        % vessel wall thickness [mm]
    PD         = pVessel.PD;        % plug flow center diameter [mm]
    profile    = pVessel.profile;   % velocity profile type [string] [parabolic0 to parabolic1: linear to quadratic]
    vMax       = pVessel.vMax;      % maximum (center) velocity [cm/s]
    vMean      = pVessel.vMean;     % mean velocity [cm/s]
    vFlow      = pVessel.vFlow;     % vessel flow [ml/min]


    % Define simulation grid
    if exist('pReal', 'var') && ~isempty(pReal)
        FOVx = pReal.FOVx; % Field of view in x-direction [mm]
        FOVy = pReal.FOVy; % Field of view in y-direction [mm]

        % first define data grid
        [pReal.xGrid, pReal.yGrid, pReal.rGrid, pReal.dx, pReal.dy, pReal.nSpin, pReal.spinDensity] = setGrid(FOVx, FOVy, size(pReal.data));
        % adjust to vessel center of mass
        tmpMask = pReal.rGrid<(ID/2);
        mag_flat = double(abs(pReal.data(tmpMask)));
        totalMag = sum(mag_flat);
        if totalMag > 0
            xCOM = sum(pReal.xGrid(tmpMask) .* mag_flat(:)) / totalMag;
            yCOM = sum(pReal.yGrid(tmpMask) .* mag_flat(:)) / totalMag;
        else
            xCOM = mean(pReal.xGrid(tmpMask));
            yCOM = mean(pReal.yGrid(tmpMask));
        end
        pReal.xGrid = pReal.xGrid - xCOM;
        pReal.yGrid = pReal.yGrid - yCOM;
        pReal.rGrid = sqrt(pReal.xGrid.^2 + pReal.yGrid.^2);

        % define simulation grid
        [xGrid, yGrid,     ~, pSim.dx, pSim.dy, nSpin, spinDensity] = setGrid(FOVx, FOVy, pSim.nSpin);
        xGrid = xGrid - xCOM;
        yGrid = yGrid - yCOM;
        rGrid = sqrt(xGrid.^2 + yGrid.^2);
    elseif isfield(pSim, 'FOVx') && isfield(pSim, 'FOVy') && ~isempty(pSim.FOVx) && ~isempty(pSim.FOVy) && pSim.FOVx ~= 0 && pSim.FOVy ~= 0
        FOVx = pSim.FOVx;
        FOVy = pSim.FOVy;
        [xGrid, yGrid, rGrid, pSim.dx, pSim.dy, nSpin, spinDensity] = setGrid(FOVx, FOVy, pSim.nSpin);
    else
        FOVx = ID+2*WT; % Field of view in x-direction [mm]
        FOVy = ID+2*WT; % Field of view in y-direction [mm]
        [xGrid, yGrid, rGrid, pSim.dx, pSim.dy, nSpin, spinDensity] = setGrid(FOVx, FOVy, pSim.nSpin);
    end
    % adjust voxel size to an odd number of simulated spins for easy definition of voxel center
    pSim.voxSizeX = ( 2*floor((pSim.voxSizeX./pSim.dx)/2) + 1 )  .*  pSim.dx;
    pSim.voxSizeY = ( 2*floor((pSim.voxSizeY./pSim.dy)/2) + 1 )  .*  pSim.dy;

    % Define masks
    % on simulation grid
    lumenMask      = rGrid<=(ID/2); % vessel lumen
    fLumen         = mean(lumenMask,[1 2]);
    wallMask       = rGrid> (ID/2) & rGrid<=(ID/2+WT); % vessel wall
    fWall          = mean(wallMask,[1 2]);
    surroundMask   = rGrid> (ID/2+WT); % static surround
    fSurround      = mean(surroundMask,[1 2]);
    lumenPlugMask  = rGrid<(PD/2) & lumenMask; % plug flow center
    fLumenPlug     = mean(lumenPlugMask,[1 2]);
    lumenLamiMask  = rGrid>=(PD/2) & lumenMask; % laminar flow region
    fLumenLami     = mean(lumenLamiMask,[1 2]);
    if exist('pReal', 'var') && ~isempty(pReal)
        % on real data grid
        pReal.lumenMask     = imresize(double(lumenMask), size(pReal.xGrid), 'box');
        pReal.wallMask      = imresize(double(wallMask), size(pReal.xGrid), 'box');
        pReal.surroundMask  = imresize(double(surroundMask), size(pReal.xGrid), 'box');
        pReal.lumenPlugMask = imresize(double(lumenPlugMask), size(pReal.xGrid), 'box');
        pReal.lumenLamiMask = imresize(double(lumenLamiMask), size(pReal.xGrid), 'box');
    end

    % Define S for each compartment
        % T_i = transverse signal magnitude of single spin in compartment i;
        % S_i = T_i * nSpin in voxel; (signal one would measure in a single-compartment voxel)
        % S_i = T_i * voxel area * spin density; (signal one would measure in a single-compartment voxel)
        % I_i = f_i * S_i; (measured signal in voxel)
    if exist('pReal', 'var') && ~isempty(pReal)
        % based on real data
        IfReal = sum(abs(pReal.data(:)).*pReal.lumenMask(:)   );
        IwReal = sum(abs(pReal.data(:)).*pReal.wallMask(:)    );
        IsReal = sum(abs(pReal.data(:)).*pReal.surroundMask(:));

        if isempty(pVessel.SlumenLami)
            pVessel.SlumenLami = IfReal./fLumen;
        else
            warning('SlumenLami specified, using it instead of real data');
        end
        if isempty(pVessel.SlumenPlug)
            pVessel.SlumenPlug = IfReal./fLumen;
        else
            warning('SlumenPlug specified, using it instead of real data');
        end
        if isempty(pVessel.Swall)
            pVessel.Swall      = IwReal./fWall;
        else
            warning('Swall specified, using it instead of real data');
        end
        if isempty(pVessel.Ssurround)
            pVessel.Ssurround  = IsReal./fSurround;        
        else
            warning('Ssurround specified, using it instead of real data');
        end
    end

    % if exist('pReal', 'var') && ~isempty(pReal) && isfield(pReal, 'magProfile') && ~isempty(pReal.magProfile)
    %     % override SlumenLami and SlumenPlug
    %     IfReal = sum(abs(pReal.data(:)).*pReal.lumenMask(:)   );
    % end

    SlumenPlug = pVessel.SlumenPlug; % Plug flow signal magnitude (in an hypothetical single-compartment voxels, independent of simulation spin density)
    SlumenLami = pVessel.SlumenLami; % Lumen signal magnitude (in a hypothetical single-compartment voxel, independent of simulation spin density). Can be a length-2 vector for velocity-dependent S (inflow)
    Swall      = pVessel.Swall;     % Wall signal magnitude (in a hypothetical single-compartment voxel, independent of simulation spin density)
    Ssurround  = pVessel.Ssurround; % Static surround signal magnitude (in a hypothetical single-compartment voxel, independent of simulation spin density)



    % Define spin velocity map
    vMap = getVelMap(rGrid, ID, profile, PD); % [cm/s]
    if ~isempty(vMax) && isempty(vMean)
        vMap = scale2maxVel(vMap, vMax ); % to the desired maximum velocity
    elseif ~isempty(vMean) && isempty(vMax)
        vMap = scale2meanVel(vMap, vMean, lumenMask); % to the desired mean velocity
    else
        error('Either vMax or vMean must be specified');
    end



    % Define signal magnitude map
    magMap = zeros(size(rGrid));
    if exist('pReal', 'var') && ~isempty(pReal) && isfield(pReal, 'magProfile') && ~isempty(pReal.magProfile)
        if fLumenLami
            magMap(lumenLamiMask) = pReal.magProfile.p(rGrid(lumenLamiMask));
        end
        if fLumenPlug
            magMap(lumenPlugMask) = pReal.magProfile.p(rGrid(lumenPlugMask));
        end
        magMap(rGrid<pReal.magProfile.rPeak) = pReal.magProfile.p(pReal.magProfile.rPeak);
        magMap(lumenMask) = magMap(lumenMask)./mean(magMap(lumenMask));
        if fLumenLami
            magMap(lumenLamiMask) = magMap(lumenLamiMask).*SlumenLami./nSpin;
        end
        if fLumenPlug
            magMap(lumenPlugMask) = magMap(lumenPlugMask).*SlumenPlug./nSpin;
        end
        % magMap(lumenMask) = pReal.magProfile.p(rGrid(lumenMask));
        % magMap(rGrid<pReal.magProfile.rPeak) = pReal.magProfile.p(pReal.magProfile.rPeak);
        % magMap(lumenMask) = magMap(lumenMask)./sum(magMap(lumenMask));

        % IfReal./fLumenLami./nSpin
        %     SlumenLami    ./nSpin
    else
        if length(SlumenLami) == 1
            magMap(lumenLamiMask) = SlumenLami./nSpin; % laminar flow region of the vessel lumen (single spin magnitudes)
        else
            magMap(lumenLamiMask) = vMap(lumenLamiMask)./(2*pVessel.vMean);
            magMap(lumenLamiMask) = magMap(lumenLamiMask).*(pVessel.SlumenLami(2)-pVessel.SlumenLami(1));
            magMap(lumenLamiMask) = magMap(lumenLamiMask) + pVessel.SlumenLami(1);
            magMap(lumenLamiMask) = magMap(lumenLamiMask)./nSpin;
        end
        magMap(lumenPlugMask) = SlumenPlug./nSpin; % plug flow center of the vessel lumen (single spin magnitudes)
    end
    magMap(wallMask)      = Swall./nSpin; % wall of the vessel (single spin magnitudes)
    magMap(surroundMask)  = Ssurround./nSpin; % static surround of the vessel (single spin magnitudes)




    % Define velocity encoding parameters
    mVenc = pVenc.mVenc; % velocity encoding bandwidth [cm/s]
    dVenc = pVenc.dVenc; % velocity encoding resolution [cm/s]
    if dVenc==inf
        venc = mVenc;
        M1   = vencToM1(venc);
    else
        M1   = linspace(0, vencToM1(dVenc), round(vencToM1(dVenc)/vencToM1(mVenc))+1);
        venc = M1toVenc(M1);
    end


    if all(pSim.voxSizeX ~= pSim.FOVx) && all(pSim.voxSizeY ~= pSim.FOVy)
        % Compute bulk complex-valued signal in FOV peak voxel
        [IpeakVox,fpeakVox,IoverlapVox,foverlapVox,NoverlapVox,I,If,Is,peakMap,overlapMap,xStart,yStart] = getBulkSignalVoxelized(magMap, vMap, venc, pSim, lumenMask,[],[],noAverageFlag); %(per unit area OR per voxel, depending on pVessel.S*)
        % If = getBulkSignalVoxelized(magMap.*lumenMask    , vMap.*lumenMask    , venc, pSim, lumenMask); %(per unit area OR per voxel, depending on pVessel.S*)
        % Is = getBulkSignalVoxelized(magMap.*(1-lumenMask), vMap.*(1-lumenMask), venc, pSim, lumenMask); %(per unit area OR per voxel, depending on pVessel.S*)

    else
        % Compute bulk complex-valued signal in FOV
        I  = getBulkSignal(magMap               , vMap               , venc); %(per unit area OR per voxel, depending on pVessel.S*)
        If = getBulkSignal(magMap.*lumenMask    , vMap.*lumenMask    , venc); %(per unit area OR per voxel, depending on pVessel.S*)
        Is = getBulkSignal(magMap.*(1-lumenMask), vMap.*(1-lumenMask), venc); %(per unit area OR per voxel, depending on pVessel.S*)
    end



    %% Output results
    % grids, maps and masks
    res.xGrid       = xGrid;
    res.yGrid       = yGrid;
    res.dx          = pSim.dx;
    res.dy          = pSim.dy;
    res.rGrid       = rGrid;
    res.vMap        = vMap;
    res.magMap      = magMap;
    res.lumenMask     = lumenMask;
    res.lumenPlugMask = lumenPlugMask;
    res.lumenLamiMask = lumenLamiMask;
    res.wallMask      = wallMask;
    res.surroundMask  = surroundMask;
    % simulated parameter
    res.FOVx        = FOVx;
    res.FOVy        = FOVy;
    res.spinDensity = spinDensity;
    res.vLumenMean  = mean(vMap(lumenMask));
    res.vLumenMax   = max( vMap(lumenMask));
    % res.sLumen      = mean(magMap(lumenMask));
    % res.sWall       = mean(magMap(wallMask));
    % res.sSurround   = mean(magMap(surroundMask));
    res.fLumen      = fLumen;
    res.fLumenPlug  = fLumenPlug;
    res.fLumenLami  = fLumenLami;
    res.fWall       = fWall;
    res.fSurround   = fSurround;
    % simulated results
    res.I           = I;
    res.If          = If;
    res.Is          = Is;
    res.venc        = venc;
    res.M1          = M1;
    if exist('IpeakVox', 'var') && ~isempty(IpeakVox)
        res.IpeakVox    = IpeakVox;
        res.fpeakVox    = fpeakVox;
        res.IoverlapVox = IoverlapVox;
        res.foverlapVox = foverlapVox;
        res.NoverlapVox = NoverlapVox;
        res.voxSizeX    = permute(pSim.voxSizeX,[1 3 2]);
        res.voxSizeY    = permute(pSim.voxSizeY,[1 3 2]);
        res.xStart      = xStart;
        res.yStart      = yStart;
        res.peakMap     = peakMap;
        res.overlapMap  = overlapMap;
    end
    % input parameters
    res.pVessel     = pVessel;
    res.pVenc       = pVenc;
    res.pSim        = pSim;

    if exist('pReal', 'var') && ~isempty(pReal)
        res.pReal     = pReal;
    end
end