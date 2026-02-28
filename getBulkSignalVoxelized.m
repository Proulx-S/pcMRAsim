function [IpeakVox,fpeakVox,IoverlapVox,foverlapVox,NoverlapVox,I,If,Is,peakMap,overlapMap,xStart,yStart] = getBulkSignalVoxelized(magMap, vMap, venc, pSim, lumenMask, xStartRpl, yStartRpl, noAverageFlag)
% magMap:  single spin magnitude map
% vMap:    single spin velocity map [cm/s]
% venc:    row vector of velocity encoding values [cm/s]
% I:       row vector of bulk complex-valued signal [complex]
% noAverageFlag: boolean flag to indicate whether to average overlapping voxels

if nargin < 8
    noAverageFlag = false;
end

if ischar(pSim.voxCenterX) && contains(pSim.voxCenterX, 'rand')
    nMonteCarlo = str2double(replace(pSim.voxCenterX, 'rand', ''));
else
    nMonteCarlo = 1;
end

I           = zeros(1          ,length(venc),1                    );
If          = zeros(1          ,length(venc),1                    );
Is          = zeros(1          ,length(venc),1                    );
IpeakVox    = zeros(nMonteCarlo,length(venc),length(pSim.voxSizeX));
fpeakVox    = zeros(nMonteCarlo,length(venc),length(pSim.voxSizeX));
if noAverageFlag
    IoverlapVox = cell(nMonteCarlo,length(venc),length(pSim.voxSizeX));
    foverlapVox = cell(nMonteCarlo,1           ,length(pSim.voxSizeX));
    NoverlapVox = cell(nMonteCarlo,1           ,length(pSim.voxSizeX));
else
    IoverlapVox = zeros(nMonteCarlo,length(venc),length(pSim.voxSizeX));
    foverlapVox = zeros(nMonteCarlo,1           ,length(pSim.voxSizeX));
    NoverlapVox = zeros(nMonteCarlo,1           ,length(pSim.voxSizeX));
end
overlapMap  = cell( nMonteCarlo,1           ,length(pSim.voxSizeX));
peakMap     = cell( nMonteCarlo,1           ,length(pSim.voxSizeX));


% Loop over velocity encoding values
for iVenc = 1:length(venc)
    % Generate spin signal map
    Imap = repmat(magMap,[size(vMap,3) 1])  .*   reshape(  exp( 1i* getVelocityEncodedPhase(vMap(:),venc(iVenc)) )  , size(vMap));

    % Get reference signals
    I( 1,iVenc) = sum(  Imap            ,[1 2]);
    If(1,iVenc) = sum(  Imap(lumenMask) ,[1 2]);
    Is(1,iVenc) = sum(  Imap(~lumenMask),[1 2]);
end




% Loop over voxel sizes
xStart = nan(nMonteCarlo,1,length(pSim.voxSizeX));
yStart = nan(nMonteCarlo,1,length(pSim.voxSizeX));
for iVoxelSize = 1:length(pSim.voxSizeX)
    dVx = round(pSim.voxSizeX(iVoxelSize)./pSim.dx);
    dVy = round(pSim.voxSizeY(iVoxelSize)./pSim.dy);

    iMonteCarlo = 1; noVoxelCounter = 0;
    while iMonteCarlo <= nMonteCarlo
        
        % Initialize voxel grid position
        if exist('xStartRpl', 'var') && ~isempty(xStartRpl) && exist('yStartRpl', 'var') && ~isempty(yStartRpl)
            xStart(iMonteCarlo,1,iVoxelSize) = xStartRpl(iMonteCarlo,1,iVoxelSize);
            yStart(iMonteCarlo,1,iVoxelSize) = yStartRpl(iMonteCarlo,1,iVoxelSize);
        else
            if iMonteCarlo == 1
                % Start voxel grid so it will be centered on the vessel cross-section
                xCenter = size(magMap,2)/2+0.5;
                xStart(iMonteCarlo,1,iVoxelSize) = floor(xCenter - dVx/2);
                xStart(iMonteCarlo,1,iVoxelSize) = xStart(iMonteCarlo,1,iVoxelSize) - floor(xStart(iMonteCarlo,1,iVoxelSize)/dVx)*dVx;
                % xStart = ceil(xCenter - floor((xCenter - dVx/2)/dVx)*dVx);
                % % xStart = floor(round((xCenter - dVx/2))/dVx)*dVx+1;
                % % % xStart = floor(round((xCenter - dVx/2))/dVx);
                % % % xStart = max(xStart, 1)*dVx;

                yCenter = size(magMap,1)/2+0.5;
                yStart(iMonteCarlo,1,iVoxelSize) = floor(yCenter - dVy/2);
                yStart(iMonteCarlo,1,iVoxelSize) = yStart(iMonteCarlo,1,iVoxelSize) - floor(yStart(iMonteCarlo,1,iVoxelSize)/dVy)*dVy;
                % yStart = ceil(yCenter - floor((yCenter - dVy/2)/dVy)*dVy);
                % % yStart = floor(round((yCenter - dVy/2))/dVy)*dVy+1;
                % % % yStart = floor(round((yCenter - dVy/2))/dVy);
                % % % yStart = max(yStart, 1)*dVy;
            else
                % Start voxel grid at a random position
                xStart(iMonteCarlo,1,iVoxelSize) = randperm(dVx,1);
                yStart(iMonteCarlo,1,iVoxelSize) = randperm(dVy,1);
            end
        end

        % Create grid indices
        xVi = zeros(1,size(magMap,2)); i = 0;
        while xStart(iMonteCarlo,1,iVoxelSize)+dVx*(i+1)-1 <= size(magMap,2)
            xVi(xStart(iMonteCarlo,1,iVoxelSize)+dVx*i:xStart(iMonteCarlo,1,iVoxelSize)+dVx*(i+1)-1) = i+1;
            i = i + 1;
        end
        yVi = zeros(1,size(magMap,1)); i = 0;
        while yStart(iMonteCarlo,1,iVoxelSize)+dVy*(i+1)-1 <= size(magMap,1)
            yVi(yStart(iMonteCarlo,1,iVoxelSize)+dVy*i:yStart(iMonteCarlo,1,iVoxelSize)+dVy*(i+1)-1) = i+1;
            i = i + 1;
        end

        % Define peak voxel
        voxGridSz = [nnz(yVi)./dVy nnz(xVi)./dVx];
        peakVox = false(voxGridSz);
        peakVox(yVi(floor(length(yVi)/2)+1), xVi(floor(length(xVi)/2)+1)) = true;

        % Get voxel volume fraction
        fVox = imresize(double(lumenMask(logical(yVi), logical(xVi))), voxGridSz, 'box');
        if noAverageFlag
            foverlapVox{iMonteCarlo,1,iVoxelSize} = fVox(fVox~=0);
            NoverlapVox{iMonteCarlo,1,iVoxelSize} = fVox(:)~=0;
        else
            foverlapVox(iMonteCarlo,1,iVoxelSize) = mean(fVox(fVox   ~=0));
            NoverlapVox(iMonteCarlo,1,iVoxelSize) = sum( fVox(:)     ~=0 );
        end
        fpeakVox(   iMonteCarlo,1,iVoxelSize) =      fVox(peakVox    );

        % Mask in spin map space
        overlapMap{iMonteCarlo,1,iVoxelSize} = false(size(lumenMask));
        overlapMap{iMonteCarlo,1,iVoxelSize}(logical(yVi), logical(xVi)) = imresize(fVox~=0,[nnz(yVi) nnz(xVi)],'box');
        peakMap{iMonteCarlo,1,iVoxelSize} = false(size(lumenMask));
        peakMap{iMonteCarlo,1,iVoxelSize}(logical(yVi), logical(xVi)) = imresize(peakVox,[nnz(yVi) nnz(xVi)],'box');
        
        if any(lumenMask(~overlapMap{iMonteCarlo,1,iVoxelSize}))
            error(['overlapping voxels do not fully cover the vessel cross-section' newline 'increase FOV to 3*ID']);
        end

        for iVenc = 1:length(venc)
            % Generate spin signal map
            Imap = repmat(magMap,[size(vMap,3) 1])  .*   reshape(  exp( 1i* getVelocityEncodedPhase(vMap(:),venc(iVenc)) )  , size(vMap));
        
            % Resample spin signal map to voxel grid
            Ivox = imresize(Imap(logical(yVi), logical(xVi), :), voxGridSz, 'box');
            Ivox = Ivox .* (dVy.*dVx); % scale back to the sum of spins            

            % Get signal from peak voxel
            IpeakVox(iMonteCarlo,iVenc,iVoxelSize) = Ivox(peakVox);

            % Get signal from overlapping voxels
            if noAverageFlag
                IoverlapVox{iMonteCarlo,iVenc,iVoxelSize} = Ivox(fVox   ~=0);
            else
                IoverlapVox(iMonteCarlo,iVenc,iVoxelSize) = sum( Ivox(fVox   ~=0));
            end
        end

        iMonteCarlo = iMonteCarlo + 1;
    end
end
