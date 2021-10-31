function [dynamicRec,staticRec,err,objectTotal,maskDynamic,maskStatic] = ...
    insitu_full_simulator(powerScale,structure,noiseLevel,propDist,binning,Ni,plotFreq,whichAlgorithm,whichDenoising)
%
% This function takes acoustic COMSOL simulations to create a dynamic
% grating. It then illuminates this grating with two D shaped probes.
% Finally, the object is recontructed using in-situ CDI.
%
% Outputs:
%
%       dynamicRec                  The reconstructed dynamic region
%
%       staticRec                   The reconstructed static region
%
%       err                         Reconstruction error as a function of
%                                   iteration number
%
%       objectTotal                 Ground truth object used to generate
%                                   the simulation
%
%       maskDynamic                 Mask indicaing which region is dynamic
%
%       maskStatic                  Mask indicating which region is static
%
%
% Inputs:
%
%       powerScale                  Ratio of the total power in the static
%                                   region to the total power in the
%                                   dynamic region
%
%       structure                   string that describes what type of
%                                   structure to put on the static region.
%                                   Options are 'flat', 'speckle',
%                                   'cameraman', 'vertGrating',' and
%                                   'horGrating'
%
%       noiseLevel                  which level of readout noise to use.
%                                   Options are 0, 1, 2, 3, & 4
%
%       propDist                    Distance to propagate the probe from
%                                   the imaging plane to the sample
%
%       binning                     Binning factor for the detector.
%
%       Ni                          number of iterations to run the
%                                   reconstruction algorithm
%
%       plotFreq                    frequency to update the plot of the
%                                   reconstruction progress. Set to Inf to
%                                   never plot
%
%       whichAlgorithm              which reconstruction algorithm to use
%                                   for the reconstruction. Options are
%                                   'inSitu', and 'inSituRAAR'
%
%       whichDenoising              which de-noising technique to use.
%                                   Options are 'BM4D', 'Thresholding', and
%                                   'None'

%% default values for inputs
if isempty(propDist)
    propDist = 0;
end
if isempty(Ni)
    Ni = 10000;
end
if isempty(plotFreq)
    plotFreq = inf;
end
if isempty(whichAlgorithm)
    whichAlgorithm = 'inSituRAAR';
end
if isempty(whichDenoising)
    whichDenoising = 'None';
end

%% Add required paths

addpath('Ptychography Data Generation Code');
addpath('BM4D_v3p2 Matlab Source Code');
addpath('InSituCDI with Noise');
%% Generate the object and probe
[objectTotal,probeTotal,maskStatic,maskDynamic] = generateGrating_JoshSim();


%% Change the static structure

% structure = 'flat'; % Options are 'flat' 'speckle' 'vertGrating' 'horGrating' and 'cameraman'

switch structure
    case 'flat'
        temp = ones(size(probeTotal));
    case 'speckle'
        temp = probeTotal.*rand(size(probeTotal));
    case 'cameraman'
        temp = imread('cameraman.tif');
        temp = double(temp);
        temp = temp/255;
        [Ny,Nx] = size(probeTotal);
        [Nyc,Nxc] = size(temp);
        temp = repmat(temp,[Ny/Nyc,Nx/Nxc]);
    case 'vertGrating'
        % Grating lines are vertical
        [Ny,Nx] = size(probeTotal);
        temp = ones(Ny,Nx);
        lambda = 28.9e-9;
        D2Det = 4E-2;
        detPix = 13.5e-6;
        Npix = 2048;

        pix = lambda*D2Det/(Npix*detPix);
        
        period = 5e-6;
        dutyCycle = 0.25;
        
        [xx,~] = meshgrid((-Nx/2:Nx/2-1)*pix,(-Ny/2:Ny/2-1)*pix);

        maskLine = mod(xx,period)>0 & mod(xx,period)<period*dutyCycle;
        maskSub = 1-maskLine;
        
        temp(maskLine==1) = 0.5;
        temp(maskSub==1) = 0.25;
    case 'horGrating'
        % Grating lines are horizontal
        [Ny,Nx] = size(probeTotal);
        temp = ones(Ny,Nx);
        lambda = 28.9e-9;
        D2Det = 4E-2;
        detPix = 13.5e-6;
        Npix = 2048;

        pix = lambda*D2Det/(Npix*detPix);
        
        period = 5e-6;
        dutyCycle = 0.25;
        
        [~,yy] = meshgrid((-Nx/2:Nx/2-1)*pix,(-Ny/2:Ny/2-1)*pix);

        maskLine = mod(yy,period)>0 & mod(yy,period)<period*dutyCycle;
        maskSub = 1-maskLine;
        
        temp(maskLine==1) = 0.5;
        temp(maskSub==1) = 0.25;
end

probeTotal(maskStatic==1) = probeTotal(maskStatic==1).*temp(maskStatic==1);

%% Change the power scaling

% powerScale = 1;

powerStatic = sum(sum(abs(probeTotal.*objectTotal(:,:,1).*maskStatic).^2));
[~,~,Nf] = size(objectTotal);

powerDynamic = zeros(1,Nf);
for ff = 1:Nf
    powerDynamic(ff) = sum(sum(abs(probeTotal.*objectTotal(:,:,ff).*maskDynamic).^2));
end

probeTotal(maskStatic==1) = probeTotal(maskStatic==1).*mean(powerDynamic)/powerStatic;


% powerScale = 10;
probeTotal(maskStatic==1) = powerScale.*probeTotal(maskStatic==1);

%% initialization
[camSetting, probe, object, scanGrid, data] = initializeInputs_v01;


%% modify simulation parameters

% =================== begin probe ===================

probe.prb0 = probeTotal;
[Ny,~] = size(probe.prb0);
% probe.sigma = 3;

% % adjust photon flux
probe.nPh_s = 1E9;
probe.lambda = 28.9e-9;
probe.D2D = 4E-2;

% =================== begin scanGrid ===================

scanGrid.scanStepSize = 12; % This means R/this number: size of ptych step

% =================== begin camSetting ===================

if noiseLevel == 0
    camSetting.pnsFlag = 0; %Poisson noise
    camSetting.wgnFlag = 0; % White gaussian noise
    camSetting.readRate = 0; % Readout rate noise: 0-4 depending on which rate
else
    camSetting.pnsFlag = 1;
    camSetting.wgnFlag = 0;
    camSetting.readRate = noiseLevel;
end

camSetting.N = Ny;
camSetting.Dt_exp = 15; %seconds?
camSetting.numAcc = 3;
camSetting.bin = binning;


%% Propagate probe slightly

% propDist = 1e-6;

camSetting.dxd = 13.5e-6;
probe.dxs = probe.D2D*probe.lambda/camSetting.N/camSetting.dxd;

if propDist ~= 0
    
    probe.prb0 = propagator_v02(probe.prb0,probe.dxs,probe.lambda,propDist,[],[]);
    
    blurWidth = 10;
    filterSize = 10;
    kern = fspecial('Gaussian',filterSize,blurWidth);
    maskStatic = filter2(kern,maskStatic,'same');
    maskStatic = maskStatic>0;

    maskDynamic = filter2(kern,maskDynamic,'same');
    maskDynamic = maskDynamic>0;
end


%% Generate a diffraction pattern for each time frame

[Ny,Nx,Nf] = size(objectTotal);
dataFull = zeros(Ny/camSetting.bin,Nx/camSetting.bin,Nf);
dataFull_clean = zeros(Ny/camSetting.bin,Nx/camSetting.bin,Nf);

for ff = 1:Nf

    % =================== begin object ===================

    object.obj0 = objectTotal(:,:,ff);

    % generate diffraction data
    [camSetting, probe, object, scanGrid, data] = ...
        GenPtyScanFFTShift_v01(camSetting, probe,object, scanGrid, data);
    
    dataFull(:,:,ff) = data.I;
    dataFull_clean(:,:,ff) = data.I_avg;
end

%% change the scaling of the data
% vacuum permittivity
eps0 = 8.8542e-12; % A^2 * s^4/kg/m^3 OR F/m

% vacuum permeability
mu0 = 1.2566e-6; % N/A^2

% speed of light
c = sqrt(1/eps0/mu0); % m/s


% undoing the physical scaing fro  m the realistic simulation to match
% propagator in the reconstruction code.

% Factors from the probe normalization (is squared later because these were
% imposed in amplitude)
probeNorm = sqrt( 2*camSetting.hw/(c*eps0) * probe.nPh_s/sum( abs(probeTotal(:)).^2 * probe.dxs * probe.dys ) );
% Factors from the ESW propagator
eswNorm = (c*eps0/2)/probe.lambda^2/probe.D2D^2*probe.dxs^2*probe.dys^2*(camSetting.N^2);
% Factors from the camera exposure
camNorm = camSetting.QE*camSetting.Dt_exp*camSetting.dxd^2/camSetting.hw;

dataFull = dataFull/(probeNorm^2)/eswNorm/camNorm*camSetting.bin^2; 
dataFull_clean = dataFull_clean/(probeNorm^2)/eswNorm/camNorm*camSetting.bin^2; 
% dataFull = dataFull/eswNorm;

%% Workaround noise simulation

% I want to apply noise levels directly to the measured diffraction pattern
% so that I can better understand the mechanism at play here.
workaround = 0;
if workaround == 1
        % convert to units of electrons (from counts)
    dataFull_clean_electron = dataFull_clean*camSetting.dynamicRange/(2^camSetting.bitDepth)*camSetting.Dt_exp;
    % Apply Poisson noise
    dataFull_electron = poissrnd(dataFull_clean_electron*100);
    dataFull = dataFull_electron /camSetting.dynamicRange*(2^camSetting.bitDepth)/camSetting.Dt_exp/100;
end

%% diffraction pattern denoising
switch whichDenoising
    case 'BM4D'
        %% BM4D denoising
        
        % figure out the noise statistics, i.e., the standard deviation of the Gaussian noise
        sigma = std(dataFull_clean(:) - dataFull(:));

        % parameters for BM4D
        distribution      = 'Gauss'; % noise distribution
                                     %  'Gauss' --> Gaussian distribution
                                     %  'Rice ' --> Rician Distribution
        profile           = 'mp';    % BM4D parameter profile
                                     %  'lc' --> low complexity
                                     %  'np' --> normal profile
                                     %  'mp' --> modified profile
                                     % The modified profile is default in BM4D. For 
                                     % details refer to the 2013 TIP paper.
        do_wiener         = 1;       % Wiener filtering
                                     %  1 --> enable Wiener filtering
                                     %  0 --> disable Wiener filtering
        verbose           = 1;       % verbose mode


        % call BM4D code
        [dataFull_denoise, ~] = ...
            bm4d(dataFull, distribution, sigma, profile, do_wiener, verbose);

        % % print information
        % disp(['==============\n input sigma = ', num2str(sigma),'\n output sigma = ', num2str(sigma_est(1,1,1))]);

        %% save denoised data (BM4D takes 24min for a 2048x2048x9 matrix!!!)
        fp = 'InSituCDI with Noise/';
        fn = ['denoisedData_noiseLevel',num2str(noiseLevel),'_',structure,'StaticStruct_powerScale',num2str(powerScale),'.mat'];
        save(fullfile(fp,fn), 'dataFull_denoise');
        
    case 'Thresholding'
        threshold = 0.2;
        dataFull_denoise = dataFull .* (dataFull > threshold);
        
    case 'None'
        dataFull_denoise = dataFull;
end

% %% load BM4D denoising result
% load('denoisedData_noiseLevel1_flatStat_powerScale1.mat', 'dataFull_denoise')

%% Reconstruct the data using whichAlgorithm ('inSitu' or 'inSituRAAR')

% in situ parameters
gamma = 0.5;

% Select a sub region based on the binning
objectTotal = objectTotal(Ny/2-Ny/2/camSetting.bin+1:Ny/2+Ny/2/camSetting.bin,Nx/2-Nx/2/camSetting.bin+1:Nx/2+Nx/2/camSetting.bin,:);
probeTotal = probeTotal(Ny/2-Ny/2/camSetting.bin+1:Ny/2+Ny/2/camSetting.bin,Nx/2-Nx/2/camSetting.bin+1:Nx/2+Nx/2/camSetting.bin);
maskStatic = maskStatic(Ny/2-Ny/2/camSetting.bin+1:Ny/2+Ny/2/camSetting.bin,Nx/2-Nx/2/camSetting.bin+1:Nx/2+Nx/2/camSetting.bin);
maskDynamic = maskDynamic(Ny/2-Ny/2/camSetting.bin+1:Ny/2+Ny/2/camSetting.bin,Nx/2-Nx/2/camSetting.bin+1:Nx/2+Nx/2/camSetting.bin);
probe.sphs0 = probe.sphs0(Ny/2-Ny/2/camSetting.bin+1:Ny/2+Ny/2/camSetting.bin,Nx/2-Nx/2/camSetting.bin+1:Nx/2+Nx/2/camSetting.bin);

static = objectTotal(:,:,1).*maskStatic;
dynamic = objectTotal.*maskDynamic;



%%

switch whichAlgorithm
    case 'inSitu'
        [dynamicRec,staticRec,err] = inSituRecon_GPU_active(sqrt(dataFull_denoise),maskStatic,maskDynamic,static,probeTotal,[],Ni,probe.sphs0,gamma,plotFreq);
        
    case 'inSituRAAR'
        % RAAR parameters
        beta = 1;    % step size
        eta1 = 0.8;    % modulus reflection depth
        eta2 = -1;    % overlap reflection depth
        
        % Run inSituRAAR
        [dynamicRec,staticRec,err] = inSitu_RAAR_Recon_GPU_active_v04(sqrt(dataFull_denoise),...
            maskStatic,maskDynamic,static,probeTotal,[],Ni,probe.sphs0,gamma,beta,eta1,eta2,plotFreq);
        
    case 'inSituUCLA'
        % add the path to the UCLA code
        addpath UCLA_code/codes
        if plotFreq==inf
            showprogress = 0;
        else
            showprogress = 1;
        end
        option=5;%custom cropping for plotting object
        
        [Recon, err, ~] = isCDI(sqrt(dataFull_denoise), probeTotal, maskStatic, Ni, showprogress, option);
        staticRec = Recon.*maskStatic;
        dynamicRec = Recon.*maskDynamic;
end
end
