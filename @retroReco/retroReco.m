classdef retroReco
    
    % Reconstruction and movie class for retrospective app
    
    properties
        
        movieExp                        % movie for movie export
        movieApp                        % movie for viewing in the app
        senseMap                        % sense map
        sos = true                      % Sum of squares reco true or false
        rescaleSlope                    % Dicom info RescaleSlope for image scaling
        rescaleIntercept                % Dicom info RescaleIntercept for image scaling
        multiSliceFlag = false          % multi-slice true or false
        multiDynamicFlag = false        % mutli-dynamic true or false
        
    end
    
    
    
    
    
    methods
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroReco
            
        end % retroReco
        
        
        
        % ---------------------------------------------------------------------------------
        % 2D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco2D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % Perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco2D_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            else
                
                % Perform CS reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                drawnow;
                cs_reco2D_mat_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mat_mc
                
                kSpaceIn = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;
                
                % kspace = {coil}[frames, x, y, slice, dynamic]
                %                   1     2  3    4       5
                dimf = size(kSpaceIn{1},1);
                dimz = size(kSpaceIn{1},4);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimf == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data x,y,frames,slices,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,1,4,5]);
                end
                
                % kspace data x,y,frames,slices,dynamics,coils
                kSpace = zeros([size(kSpaceIn{1}),nrCoils]);
                for i = 1:nrCoils
                    kSpace(:,:,:,:,:,i) = kSpaceIn{i};
                end
                
                % averages data x,y,frames,slices,dynamics
                averages = permute(averages,[2,3,1,4,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimx = 2^nextpow2(size(kSpace,1));
                dimy = 2^nextpow2(size(kSpace,2));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kSpace,3);
                
                % pre-allocate memory for image_out
                imageOut = zeros(dimf,dimx,dimy,dimz,dimd);
                
                % slice and dynamic loop
                for slice = 1:dimz
                    
                    for dynamic = 1:dimd
                        
                        % kspace of slice and dynamic
                        kData = squeeze(kSpace(:,:,:,slice,dynamic,:));
                        
                        % zero padding
                        padSizex = round((dimx - size(kData,1))/2);
                        padSizey = round((dimy - size(kData,2))/2);
                        kData = padarray(kData,[padSizex,padSizey,0],'both');
                        kData = kData(1:dimx,1:dimy,:,:);
                        
                        % size of the data
                        [nx,ny,~,nrCoils]=size(kData);
                        
                        % normalize the data in the range of approx 0 - 1 for better numerical stability
                        kData = kData/max(abs(kData(:)));
                        
                        % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                        mask = squeeze(averages(:,:,:,slice,dynamic));
                        mask = padarray(mask,[padSizex,padSizey,0],'both');
                        mask = mask(1:dimx,1:dimy,:);
                        mask = mask./mask;
                        mask(isnan(mask)) = 1;
                        mask = logical(mask);
                        
                        % coil sensitivity map
                        b1 = ones(nx,ny,nrCoils);
                        
                        % data
                        param.y = kData;
                        
                        % reconstruction design matrix
                        param.E = Emat_yxt(mask,b1);
                        
                        % Total variation (TV) constraint in the temporal domain
                        % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                        % TV only in the time domain
                        param.TV = TVOP;
                        param.TVWeight = lambdaTV/8;
                        
                        % number of iterations, 2 x 10 iterations
                        param.nite = 10;
                        param.nouter = 2;
                        param.totaliterations = dimz * dimd * param.nouter * param.nite;
                        
                        % linear reconstruction
                        kData1 = randn(size(kData))/2000 + kData;  % add a little bit of randomness, such that linear reco is not exactly right
                        reconDft = param.E'*kData1;
                        
                        % iterative reconstruction
                        reconCS = reconDft;
                        for n = 1:param.nouter
                            [reconCS,param.iteration] = CSL1NlCg(app,reconCS,param);
                        end
                        
                        % rearrange to correct orientation: frames, x, y
                        imageTmp = flip(permute(squeeze(abs(reconCS)),[3, 1, 2]),2);

                        % output reconstructed image
                        imageOut(:,:,:,slice,dynamic) = imageTmp;
                        
                    end
                    
                end
                
                % correct back to 1 frame reconstruction
                if dimf == 1
                    imageOut = imageOut(1,:,:,:,:);
                end
           
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mc
                
                kSpaceIn = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxy = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                SOS = objReco.sos;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in CINE dimension
                % TVxy = total variation in xy-dimension regularization
                % TVd = total variation in dynamic dimension
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[CINE, x, y, slice, dynamic]
                %                     1    2  3    4       5
                dimf = size(kSpaceIn{1},1);
                dimx = 2^nextpow2(size(kSpaceIn{1},2));
                dimy = 2^nextpow2(size(kSpaceIn{1},3));
                dimz = size(kSpaceIn{1},4);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize k-space to next power of 2
                for i = 1:nrCoils
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy)],kSpaceIn{i});
                end
                
                % put all data in a normal matrix
                kSpace = zeros(dimf,dimx,dimy,dimz,dimd);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimf
                            for i = 1:nrCoils
                                kSpace(j,:,:,k,l,i) = kSpaceIn{i}(j,:,:,k,l);
                            end
                        end
                    end
                end
                
                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15
                
                kspacePics = permute(kSpace,[7,3,2,6,8,9,10,11,12,13,1,5,14,4]);
                
                % wavelet in spatial dimensions 2^1+2^2=6
                % total variation in spatial dimensions 2^1+2^2=6
                % total variation in cine dimension 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048
                
                if ESPIRiT && nrCoils>1
                    
                    % ESPIRiT reconstruction
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    
                    % Calculate coil sensitivity maps with ecalib bart function
                    kspace_pics_sum = sum(kspacePics,[11,12]);
                    sensitivities = bart(app,'ecalib -S -I -a', kspace_pics_sum);      % ecalib with softsense
                    
                    picsCommand = 'pics ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picsCommand = [picsCommand,' -RT:6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspacePics,sensitivities);
                    
                    % Sum of squares reconstruction
                    if SOS
                        imageReg = bart(app,'rss 16', imageReg);
                    end
                    imageReg = abs(imageReg);
                    
                else
                    
                    % Reconstruction without sensitivity correction, including coil scaling
                    sensitivities = ones(1,dimy,dimx,nrCoils,1,1,1,1,1,1,1,1,1,dimz);
                    
                    % regular reconstruction
                    picsCommand = 'pics ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picsCommand = [picsCommand,' -RT:6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        picsCommand = [picsCommand,' -b20 -RL:6:6:',num2str(LR)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspacePics,sensitivities);

                    % total-generalized-variation denoising. Seems not really useful
                    % imageReg2 = bart(app,'tgv 0.001 1030',imageReg);
                    % imageReg = imageReg2(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,1);
                    
                    % Take absolute values
                    imageReg = abs(imageReg);
                    
                end
                
                %                              y    x  frames dynamic slices
                imageReg = reshape(imageReg,[dimy,dimx,dimf,dimd,dimz]);

                % rearrange to correct orientation: frames, x, y, slices, dynamics
                imageOut = flip(permute(imageReg,[3,2,1,5,4]),3);

                % sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[3,2,14,4,5,1,6,7,8,9,10,11,12,13]),2);
                
                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                
            end
            
        end
        
                
        
        % ---------------------------------------------------------------------------------
        % 3D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco3D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the 3D data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco3D_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            else
                
                % perform cs reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                cs_reco3D_mat_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mat_mc
                
                kSpaceIn = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;
     
                % kdata_in = {coil}[frames, x, y, z, dynamics]
                %                      1    2  3  4     5
                dimf = size(kSpaceIn{1},1);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimf == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data x,y,z,frames,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,4,1,5]);
                end

                % kspace data x,y,z,frames,dynamics,coils
                kSpace = zeros([size(kSpaceIn{1}),nrCoils]);
                for i = 1:nrCoils
                    kSpace(:,:,:,:,:,i) = kSpaceIn{i};
                end
                
                % averages data x,y,z,frames,dynamics
                averages = permute(averages,[2,3,4,1,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimx = 2^nextpow2(size(kSpace,1));
                dimy = 2^nextpow2(size(kSpace,2));
                dimz = 2^nextpow2(size(kSpace,3));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kSpace,4);
                
                % pre-allocate memory for image_out
                imageOut = zeros(dimf,dimx,dimy,dimz,dimd);
                
                for dynamic = 1:dimd
                    
                    % kspace of slice and dynamic
                    kData = squeeze(kSpace(:,:,:,:,dynamic,:));
                    
                    % zero padding
                    padSizex = round((dimx - size(kData,1))/2);
                    padSizey = round((dimy - size(kData,2))/2);
                    padSizez = round((dimz - size(kData,3))/2);
                    kData = padarray(kData,[padSizex,padSizey,padSizez,0],'both');
                    kData = kData(1:dimx,1:dimy,1:dimz,:,:);
                    
                    % size of the data, z,y,x,frames,coils
                    [nx,ny,nz,~,nrCoils] = size(kData);
                    
                    % normalize the data in the range of approx 0 - 1 for better numerical stability
                    kData = kData/max(abs(kData(:)));
                    
                    % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                    mask = squeeze(averages(:,:,:,:,dynamic));
                    mask = padarray(mask,[padSizex,padSizey,padSizez,0],'both');
                    mask = mask(1:dimx,1:dimy,1:dimz,:);
                    mask = mask./mask;
                    mask(isnan(mask)) = 1;
                    mask = logical(mask);
                    
                    % coil sensitivity map
                    b1 = ones(nx,ny,nz,nrCoils);
                    
                    % data
                    param.y = kData;
                    
                    % reconstruction design matrix
                    param.E = Emat_zyxt(mask,b1);
                    
                    % Total variation (TV) constraint in the temporal domain
                    % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                    % TV only in the time domain
                    param.TV = TVOP3D;
                    param.TVWeight = lambdaTV/8;
                    
                    % number of iterations, 2 x 10 iterations
                    param.nite = 10;
                    param.nouter = 2;
                    param.totaliterations = dimd * param.nouter * param.nite;
                    
                    % linear reconstruction
                    kData1 = randn(size(kData))/2000 + kData;  % add a little bit of randomness, such that linear reco is not exactly right
                    reconDft = param.E'*kData1;
                    
                    % iterative reconstruction
                    reconCs=reconDft;
                    for n=1:param.nouter
                        [reconCs,param.iteration] = CSL1NlCg(app,reconCs,param);
                    end
                    
                    % rearrange to correct orientation: frames, x, y, z
                    imageTmp = flip(permute(abs(reconCs),[4, 1, 2, 3]),2);

                    % output reconstructed image
                    imageOut(:,:,:,:,dynamic) = imageTmp;
                    
                end
                
                % correct back to 1 frame reconstruction
                if dimf == 1
                    imageOut = imageOut(1,:,:,:,:);
                end
                
                % shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mc
                
                kSpaceIn = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in time regularization
                % TVxyz = total variation in xyz-dimension regularization
                % VNSA = variable number of signal averages correction true or false
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[frames, x, y, z, dynamics]
                %                     1      2  3  4     5
                dimf = size(kSpaceIn{1},1);
                dimx = 2^nextpow2(size(kSpaceIn{1},2));
                dimy = 2^nextpow2(size(kSpaceIn{1},3));
                dimz = 2^nextpow2(size(kSpaceIn{1},4));
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize to next power of 2
                for i = 1:nrCoils
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy),' 3 ',num2str(dimz)],kSpaceIn{i});
                end
                
                % kspace suitable for bart
                kSpace = zeros(dimf,dimx,dimy,dimz,dimd,nrCoils);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimf
                            for i = 1:nrCoils
                                kSpace(j,:,:,k,l,i) = kSpaceIn{i}(j,:,:,k,l);
                            end
                        end
                    end
                end
                
                % Bart dimensions
                % 	READ_DIM,       1   z
                % 	PHS1_DIM,       2   y
                % 	PHS2_DIM,       3   x
                % 	COIL_DIM,       4   coils
                % 	MAPS_DIM,       5   sense maps
                % 	TE_DIM,         6
                % 	COEFF_DIM,      7
                % 	COEFF2_DIM,     8
                % 	ITER_DIM,       9
                % 	CSHIFT_DIM,     10
                % 	TIME_DIM,       11  cardiac / respiratory CINE
                % 	TIME2_DIM,      12  dynamics
                % 	LEVEL_DIM,      13
                % 	SLICE_DIM,      14  slices
                % 	AVG_DIM,        15
                
                kspace_pics = permute(kSpace,[4,3,2,6,7,8,9,10,11,12,1,5,13,14]);
                
                % wavelet in spatial dimensions 2^0+2^1+2^2=7
                % total variation in spatial dimensions 2^0+2^1+2^2=7
                % total variation in time 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048
                
                if ESPIRiT && nrCoils>1
                    
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    kspace_pics_sum = sum(kspace_pics,[11,12]);
                    sensitivities = bart(app,'ecalib -I -S -a', kspace_pics_sum);
                    
                    app.ProgressGauge.Value = 25;
                    drawnow;
                    
                    picsCommand = 'pics -S ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picsCommand = [picsCommand,' -RT:7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspace_pics,sensitivities);
                    
                    app.ProgressGauge.Value = 95;
                    drawnow;
                    imageReg = bart(app,'rss 16', imageReg);
                    imageReg = abs(imageReg);
                    
                else
                    
                    % Reconstruction without sensitivity correction
                    sensitivities = ones(dimz,dimy,dimx,nrCoils,1,1,1,1,1,1,1,1,1,1);
                    
                    picsCommand = 'pics -S ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picsCommand = [picsCommand,' -RT:7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        picsCommand = [picsCommand,' -RL:7:7:',num2str(LR)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    
                    imageReg = abs(bart(app,picsCommand,kspace_pics,sensitivities));
                    
                end
                
                % rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(imageReg,[dimz,dimy,dimx,dimf,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,3,2,1,5]),3),4);
                imageOut = circshift(imageOut,1,4);
                
                % sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);
                senseMap1 = circshift(senseMap1,1,4);
                
                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));
                
                % shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift2,4);
                
            end
            
        end % reco3D
        
                
        
        % ---------------------------------------------------------------------------------
        % 2D radial reconstruction with the Bart toolbox /// NOT FUNCTIONAL ///
        % ---------------------------------------------------------------------------------
        function objReco = reco2Dradial(objReco, objData, objKspace, app)
            
            app.TextMessage('Reconstructing 2D radial data with the BART toolbox ...');
            app.ProgressGauge.Value = 0;
            drawnow;
            
            kspace_in = objKspace.kSpace;
            Wavelet = app.WVxyzEditField.Value;
            TVxy = app.TVxyzEditField.Value;
            TVt = app.TVcineEditField.Value;
            nc = objData.nr_coils;
            dimx = objData.dimx;
            dimy = objData.dimy;
            ksteps = objData.nrKsteps; %#ok<NASGU> 
            
            % app = matlab app
            % kspace_in = sorted k-space
            % nc = number of RF receiver coils
            % averages_in = number of averages per k-line
            % Wavelet = wavelet L1-norm regularization factor
            % TVt = total variation in time regularization
            % TVxy = total variation in xy-dimension regularization
            
            % radial k-space trajectory for Bart toolbox
            radial = 360;           % 360 degrees
            u_spokes = true;        % unique spokes true/false
            half_spoke = false;     %#ok<NASGU> % half trajectory true/false
            
            for n = 1:dimy
                
                % caluculate angle
                if radial == 180
                    a(n) = (n-1)*180/dimy; %#ok<*AGROW> 
                elseif ~u_spokes
                    a(n) = (n-1)*360/dimy;
                elseif (n-1) < dimy/2
                    a(n) = (n-1)*360/dimy;
                else
                    a(n) = (n-1)*360/dimy + 180/dimy;
                end
                
                % calculate x and y values
                for m = 1:dimx
                    r(m,n) = (dimx-1)/2 -(m-1);
                    x = r(m,n)*-cos(a(n)*(pi/180));
                    y = r(m,n)*sin(a(n)*(pi/180));
                    traj(1,m,n) = x;
                    traj(2,m,n) = y;
                    traj(3,m,n) = 0;
                end
                
            end
            
            % bart (kz,ky,kx,channels,multiple sense, ..., ..., cardiac frames at pos 11)
            % sense maps for 2D data:  (kz,ky,kx,....)
            for i = 1:nc
                kspace(:,:,:,:,i) = kspace_in{i};
            end
            
            %figure(1)
            %imshow(real(squeeze(kspace(1,:,:))),[-2*pi,2*pi]);
            %imshow(abs(squeeze(kspace(1,:,:))),[0 10]);
            
            % rearrange for the bart toolbox
            kspace_pics = permute(kspace,[4,3,2,5,6,7,8,9,10,11,1]);
            
            % create trajectory
            %bartcommand = ['traj -r -y',num2str(ksteps),' -x',num2str(dimx),' -q0:0:0'];
            %traj = bart(bartcommand);
            
            % sensitivity map
            kspace_pics_sum = sum(kspace_pics,11);
            lowres_img = bart(app,'nufft -i -l6 -d16:16:1 -t', traj, kspace_pics_sum);
            lowres_ksp = bart(app,'fft -u 7', lowres_img);
            
            highres_img = abs(bart(app,'nufft -i -t -l6 -d128:128:1 -t', traj, kspace_pics_sum)); %#ok<NASGU> 
            figure(1)
            imshow(abs(lowres_img),[0,1.5*max(abs(lowres_img(:)))]);
            
            % zeropad to full size
            bartcommand = ['resize -c 0 ',num2str(dimx),' 1 ',num2str(dimx)];
            ksp_zerop = bart(app, bartcommand, lowres_ksp);
            
            % calculate sensitivity map
            sensitivities = bart(app, 'ecalib -m1', ksp_zerop); %#ok<NASGU> 
            
            sensitivities = ones(dimx,dimx,1,1,1,1,1,1,1,1,1,1,1);
            
            picscommand = ['pics -S -u1 -RW:6:0:',num2str(Wavelet),' -RT:6:0:',num2str(TVxy),' -RT:1024:0:',num2str(TVt),' -t'];
            image_reg = bart(app,picscommand,traj,kspace_pics,sensitivities); %#ok<NASGU> 
            
            image_reg = bart(app,'nufft -i -l1 -d128:128:1 -t', traj, kspace_pics);
            
            % rearrange to correct orientation: frames,x, y
            image_reg = abs(squeeze(image_reg(:,:,1,1,1,1,1,1,1,1,:)));
            image_out = flip(permute(image_reg,[3, 2, 1]),3);
            sensemap = flip(permute(squeeze(abs(sensitivities)),[2, 1, 3, 4]),2);
            
            % normalize sense map to reasonable value range
            sensemap = sensemap*1024/max(sensemap(:));
            
            % shift image in phase-encoding direction if needed
            objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
            objReco.senseMap = circshift(sensemap,objData.pixelshift1,2);
            
            app.ProgressGauge.Value = 100;
            drawnow;
            
        end % reco2Dradial
        
        
         
        
        % ---------------------------------------------------------------------------------
        % 3D UTE reconstruction with the Bart toolbox
        % ---------------------------------------------------------------------------------
        function objReco = reco3Dute(objReco, objData, objKspace, app)

            if app.bartDetectedFlag

                app.TextMessage('Reconstructing 3D UTE data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;

                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                nCoils = objData.nr_coils;

                % Wavelet = wavelet L1-norm regularization factor
                % TVxyz = total variation in xyz-dimension regularization
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVd = total variation in dynamic dimension regularization
                % nc = number of RF receiver coils
                
                dimx = size(objKspace.kSpace{1},4);
                dimy = dimx;
                dimz = dimx;
                dimf = size(objKspace.kSpace{1},1);
                dimd = size(objKspace.kSpace{1},5);
               
                kSpace = zeros(size(objKspace.kSpace{1}));
                for i = 1:nCoils
                    kSpace(:,:,:,:,:,i) = objKspace.kSpace{i};
                end
                traj = objKspace.kSpaceTraj;
                
                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15
             
                % Rearrange for BART         1  2  3  4  5  6  7  8  9 10 11 12 13 14
                kSpacePics = permute(kSpace,[7, 4, 2, 6, 3, 8, 9,10,11,12,1, 5, 13,14]);
                
                % Rearrange for BART     1  2  3  4  5  6  7  8  9 10 11 12 13 14
                trajPics = permute(traj,[6, 4, 2, 3, 7, 8, 9,10,11,12, 1, 5,13,14]);
  
                % Coil sensitivities from sum of all frames and dynamics
                if nCoils > 1
                    kSpacePicsSum = sum(kSpacePics,[11,12]);
                    trajPicsSum = sum(trajPics,[11,12]);
                    app.TextMessage('Calculating coil sensitivity maps ...');
                    lowResImage = bart(app,'nufft -i -d64:64:64 -t', trajPicsSum, kSpacePicsSum);
                    lowResKspace = bart(app,'fft -u 7', lowResImage);
                    kSpaceZeroFilled = bart(app,['resize -c 0 ',num2str(dimz),' 1 ',num2str(dimy),' 2 ',num2str(dimx)], lowResKspace);
                    sensitivities = bart(app,'ecalib -S -t0.0005 -m1', kSpaceZeroFilled);
                else
                    sensitivities = ones(dimz,dimy,dimx,1,1,1,1,1,1,1,1,1,1,1);
                end

                % Density correction
                app.TextMessage('Calculating density correction ...');
                denseOnes = ones(size(kSpacePics));
                denseTmp = bart(app,strcat('nufft -d',num2str(dimx),':',num2str(dimx),':',num2str(dimx),' -a'),trajPics,denseOnes);
                density = bart(app,'nufft ',trajPics,denseTmp);
                density = density.^(-1/3);
                density(isnan(density)) = 0;
                density(isinf(density)) = 0;

                % PICS reconstruction
                app.TextMessage('PICS reconstruction ...');
                picsCommand = 'pics -i30 -d1 ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -RT:7:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blocksize = round(dimx/16);  % Block size
                    blocksize(blocksize<4) = 4;
                    picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -RT:1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -RT:2048:0:',num2str(TVd)];
                end
                igrid = bart(app,picsCommand,'-t',trajPics,'-p',density,kSpacePics,sensitivities);
                
                % Root sum of squares over all coils
                recoImage = bart(app,'rss 8', igrid);

                % Rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(recoImage,[dimz,dimy,dimx,dimf,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,1,2,3,5]),3),4);
                
                % Sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);
                
                % Normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));
                
                % Shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift2,4);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            else

                app.TextMessage('BART toolbox not available ...');
                app.TextMessage('3D UTE reconstruction using 3D NUFFT ...');
                app.ProgressGauge.Value = 0;
                
                dimx = size(objKspace.kSpace{1},4);
                dimy = dimx;
                dimz = dimx;
                dimf = size(objKspace.kSpace{1},1);
                dimd = size(objKspace.kSpace{1},5);
                nCoils = objData.nr_coils;
                loops = dimf*dimd*nCoils;
                app.TextMessage('Slow reconstruction ...');
                app.TextMessage(strcat('Estimated reconstruction time >',{' '},num2str(loops*2),{' '},'min ...'));
                
                traj = objKspace.kSpaceTraj;
                image = zeros(dimf,dimx,dimy,dimz,dimd,nCoils);
                sensitivities = 4095*ones(dimz,dimy,dimx,1,1,1,1,1,1,1,1,1,1,1);
              
                maxit = 10;     % 0 or 1 for gridding, higher values for conjugate gradient
                damp = 0;       % Tikhonov penalty on ||x||
                weight = [];    % data weighting (optional)
                partial = 0.5;  % Tikhobov penalty on ||imag(x))||
      
                cnt = 0;

                for dynamic = 1:dimd

                    for frame = 1:dimf

                        for coil = 1:nCoils

                            app.ProgressGauge.Value = round(100*cnt/loops);
                            drawnow;

                            om = permute(squeeze(traj(frame,:,1,:,dynamic,:)),[3,2,1]);
                            obj = nufft_3d(om,dimx);

                            data = squeeze(objKspace.kSpace{coil}(frame,:,1,:,dynamic));
                            data = permute(data,[2 1]);
                            data = data(:);
                       
                            image(frame,:,:,:,dynamic,coil) = obj.iNUFT(data,maxit,damp,weight,'phase-constraint',partial);
                       
                            cnt = cnt + 1;

                        end

                    end

                end

                % Root sum of squares over coil dimenion
                image = rssq(image,6);
    
                % Shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(image,-objData.pixelshift1,3);
                objReco.movieExp = circshift(image,-objData.pixelshift2,4);
                objReco.senseMap = circshift(sensitivities,-objData.pixelshift1,3);
                objReco.senseMap = circshift(sensitivities,-objData.pixelshift2,4);
                
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            end

        end % reco3Dute



        
        % ---------------------------------------------------------------------------------
        % Normalize movie intensity
        % ---------------------------------------------------------------------------------
        function objReco = normImages(objReco,objData)
            
            % normalize the images to 2^15 range

            objReco.rescaleIntercept = 0;                           % Dicom rescale intercept
            objReco.rescaleSlope = max(objReco.movieExp(:));        % Dicom rescale slope
   
            objReco.movieExp = round(32765*objReco.movieExp/objReco.rescaleSlope);
            objReco.movieApp = objReco.movieExp;

            if objData.PHASE_ORIENTATION == 0
                objReco.movieApp = flip(objReco.movieApp,2);
            end
            
        end % normImages
        
        

        % ---------------------------------------------------------------------------------
        % Determine whether movie has multiple slices and/or dynamics
        % ---------------------------------------------------------------------------------
        function objReco = determineMultiDimensions(objReco)
            
            % multislice
            if size(objReco.movieApp,4) > 1
                objReco.multiSliceFlag = true;
            else
                objReco.multiSliceFlag = false;
            end

            % multidynamic
            if size(objReco.movieApp,5) > 1
                objReco.multiDynamicFlag = true;
            else
                objReco.multiDynamicFlag = false;
            end

        end % determineMultiDimensions
        

        

        % ---------------------------------------------------------------------------------
        % Image reconstruction: SUR files
        % ---------------------------------------------------------------------------------
        function obj = recoSurFiles(obj, surpath, suffix, mrdfilename, rprfilename)

            % SUR file names
            surfiles = [surpath, suffix, '_00###.SUR'];

            % Link with the server
            m_Recon = actxserver('recon.Application');

            set(m_Recon,'Visible',1);
            set(m_Recon,'DisplayImages',1);

            % Filenames
            set(m_Recon,'DataFile',mrdfilename);
            set(m_Recon,'RPRFile',rprfilename);
            set(m_Recon,'ImageFile',surfiles);

            % Delete old SUR files
            scmd = ['del /Q ', surpath, '*.SUR'];
            system(scmd);

            % Do the reco
            invoke(m_Recon,'Run');

            % Wait for recon to complete
            while get(m_Recon,'StatusID')~=4
                if get(m_Recon,'StatusID')==5
                    break;
                end
                pause(0.1);
            end

            % Stop the link
            invoke(m_Recon,'Quit');

        end % recoSurFiles
        

        

        % ---------------------------------------------------------------------------------
        % L-curve calculation
        % ---------------------------------------------------------------------------------
        function obj = recoLcurve(obj,objData,objKspace,app)

            % Range of lamda's
            lambda = logspace(log10(app.Lambda1EditField.Value),log10(app.Lambda2EditField.Value),app.LambdaStepsEditField.Value);
            lambda = round(lambda,3);

            % Determine which parameters will be optimized
            par = [0 0 0 0 0];
            if app.LcurveWVxyzCheckBox.Value  == 1  par(1) = 1; end %#ok<*SEPEX>
            if app.LcurveTVxyzCheckBox.Value  == 1  par(2) = 1; end
            if app.LcurveLRxyzCheckBox.Value  == 1  par(3) = 1; end
            if app.LcurveTVcineCheckBox.Value == 1  par(4) = 1; end
            if app.LcurveTVdynCheckBox.Value  == 1  par(5) = 1; end

            % Intialize progress gauge
            totalNumberOfSteps = app.LambdaIterationsEditField.Value * sum(par) * length(lambda);
            app.ProgressGauge.Value = 0;
            cnt = 0;

            % Iterations
            iter = 0;
            while (iter<app.LambdaIterationsEditField.Value) && (app.stopLcurveFlag==false)
                iter = iter + 1;

                % Wavelet xyz
                if par(1)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| WVxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Wavelet
                        app.WVxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, WVxyz, iteration =',num2str(iter),', lambda =',num2str(lambda(i))));

                        % CS reco
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm wavelet
                        WVxyz = bart(app,'cdf97 6',Lmovie);
                        y(i) = norm(WVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation spatial domain
                if par(2)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVxyz, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in spatial dimension
                        for dim1=1:size(Lmovie,1)
                            for dim4=1:size(Lmovie,4)
                                for dim5=1:size(Lmovie,5)
                                    for dim6=1:size(Lmovie,6)
                                        [Gx,Gy] = gradient(squeeze(Lmovie(dim1,:,:,dim4,dim5,dim6)));
                                        TVxyz(dim1,:,:,dim4,dim5,dim6) = sqrt(Gx.^2 + Gy.^2);
                                    end
                                end
                            end
                        end
                        y(i) = norm(TVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Low-Rank in spatial domain
                if par(3)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| LRxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.LLRxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, LRxyz, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in spatial dimension
                        for dim1=1:size(Lmovie,1)
                            for dim4=1:size(Lmovie,4)
                                for dim5=1:size(Lmovie,5)
                                    for dim6=1:size(Lmovie,6)
                                        [Gx,Gy] = gradient(squeeze(Lmovie(dim1,:,:,dim4,dim5,dim6)));
                                        TVxyz(dim1,:,:,dim4,dim5,dim6) = sqrt(Gx.^2 + Gy.^2);
                                    end
                                end
                            end
                        end
                        y(i) = norm(TVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation CINE
                if par(4)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVcine x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVcineEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVcine, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV
                        TVcine = circshift(Lmovie,1) - Lmovie;
                        y(i) = norm(TVcine(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation dynamic dimension
                if par(5)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVdyn x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVdynEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVdyn, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in dynamic dimension
                        TVtime = circshift(Lmovie,5) - Lmovie;
                        y(i) = norm(TVtime(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0; %#ok<*AGROW>
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end


                % Analyze the L-curve and determine optimal value
                if ~app.stopLcurveFlag
                    LcurveAnalysisFnc(lambda,x,y,4);
                end


            end


            % Find the point of highest curvature
            function LcurveAnalysisFnc(lambda,x,y,par)

                % Calculate the curvature and find the maximum
                xq = log(x);
                yq = log(y);
               
                curv = zeros(length(xq),1);

                for pnt = 2:length(xq)-1
                    
                    a = sqrt( (xq(pnt+1)-xq(pnt))^2   + (yq(pnt+1)-yq(pnt))^2   );
                    b = sqrt( (xq(pnt-1)-xq(pnt))^2   + (yq(pnt-1)-yq(pnt))^2   );
                    c = sqrt( (xq(pnt-1)-xq(pnt+1))^2 + (yq(pnt-1)-yq(pnt+1))^2 );
                    s = ( a + b + c ) / 2;
                    A = sqrt(s*(s-a)*(s-b)*(s-c));
                    curv(pnt) = 4*A/(a*b*c);

                end
              
                % Set lambda to the value corresponding to maximum curvature
                [~,idx] = max(curv);

                switch par

                    case 1
                        app.WVxyzEditField.Value = lambda(idx);

                    case 2
                        app.TVxyzEditField.Value = lambda(idx);

                    case 3
                        app.LLRxyzEditField.Value = lambda(idx);

                    case 4
                        app.TVcineEditField.Value = lambda(idx);

                    case 5
                        app.TVdynEditField.Value = lambda(idx);

                end

                % Indicate the value in the L-curve and model-error curve
                hold(app.LcurveFig,"on");
                loglog(app.LcurveFig,x(idx),y(idx),'o','MarkerEdgeColor',app.orange,'MarkerFaceColor',app.brightred,'LineWidth',1.5);
                hold(app.LcurveFig,"off");

                hold(app.LcurveModelErrorFig,"on");
                loglog(app.LcurveModelErrorFig,[lambda(1) lambda(end)],[x(idx) x(idx)],'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',app.blue,'LineWidth',1.5);
                hold(app.LcurveModelErrorFig,"off");

                pause(1);

            end % LcurveAnalysisFnc

        end % recoLcurve



    end % methods
    
end % retroReco

