%% insitu Noise simulator
% This function take the obejct and probe from the generator function and
% creates diffraction patterns with realistic noise!! (or no noise if you
% want that).

function [dynamicRec,staticRec,err,objectTotal,maskDynamic,maskStatic] = insitu_noise_simulator(structure,noiseLevel,propDist,Ni,plotFreq)


%% Add required paths

addpath('Ptychography Data Generation Code');
addpath('/home/imaging/binw/Github/inSituSimulation/BM4D_v3p2 Matlab Source Code');
addpath('/home/imaging/binw/Github/inSituSimulation/InSituCDI with Noise');
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

powerScale = 1;

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
    camSetting.wgnFlag = 1;
    camSetting.readRate = noiseLevel;
end

camSetting.N = Ny;
camSetting.Dt_exp = 1; %seconds?
camSetting.numAcc = 3;

%% Propagate probe slightly

if isempty(propDist)
    propDist = 1e-6;
end

% propDist = 1e-6;

camSetting.dxd = 13.5e-6;
probe.dxs = probe.D2D*probe.lambda/camSetting.N/camSetting.dxd;

probe.prb0 = propagator_v02(probe.prb0,probe.dxs,probe.lambda,propDist,[],[]);

if propDist~=0
    blurWidth = 50;
    filterSize = 100;
    kern = fspecial('Gaussian',filterSize,blurWidth);
    maskStatic = filter2(kern,maskStatic,'same');
    maskStatic = maskStatic>0;

    maskDynamic = filter2(kern,maskDynamic,'same');
    maskDynamic = maskDynamic>0;
end


%% Generate a diffraction pattern for each time frame

[Ny,Nx,Nf] = size(objectTotal);
dataFull = zeros(Ny,Nx,Nf);
dataFull_clean = zeros(Ny,Nx,Nf);

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

dataFull = dataFull/(probeNorm^2)/eswNorm/camNorm; 
dataFull_clean = dataFull_clean/(probeNorm^2)/eswNorm/camNorm; 
% dataFull = dataFull/eswNorm;

% %% diffraction patterns skipping noise simulation
% 
% norm = sqrt(camSetting.N^2);
% for ff = 1:Nf
%     dataFull(:,:,ff) = objectTotal(:,:,ff).*probeTotal.*probe.sphs0;
%     dataFull(:,:,ff) = fftshift(fft2(ifftshift(dataFull(:,:,ff))))/norm;
%     dataFull(:,:,ff) = abs(dataFull(:,:,ff)).^2;
% end

%% BM4D denoising
% 
% % figure out the noise statistics, i.e., the standard deviation of the Gaussian noise
% sigma = std(dataFull_clean(:) - dataFull(:));
% 
% % parameters for BM4D
% distribution      = 'Gauss'; % noise distribution
%                              %  'Gauss' --> Gaussian distribution
%                              %  'Rice ' --> Rician Distribution
% profile           = 'mp';    % BM4D parameter profile
%                              %  'lc' --> low complexity
%                              %  'np' --> normal profile
%                              %  'mp' --> modified profile
%                              % The modified profile is default in BM4D. For 
%                              % details refer to the 2013 TIP paper.
% do_wiener         = 1;       % Wiener filtering
%                              %  1 --> enable Wiener filtering
%                              %  0 --> disable Wiener filtering
% verbose           = 1;       % verbose mode
% 
% 
% % call BM4D code
% [dataFull_denoise, sigma_est] = ...
%     bm4d(dataFull, distribution, sigma, profile, do_wiener, verbose);
% 
% % print information
% disp(['==============\n input sigma = ', num2str(sigma),'\n output sigma = ', num2str(sigma_est(1,1,1))]);
% 
%% save denoised data (BM4D takes 24min for a 2048x2048x9 matrix!!!)
% fp = '/home/imaging/binw/Github/inSituSimulation/InSituCDI with Noise/';
% fn = ['denoisedData_noiseLevel',num2str(noiseLevel),'_flatStat_powerScale',num2str(powerScale),'.mat'];
% save([fp fn], 'dataFull_denoise');

%% load BM4D denoising result
load('denoisedData_noiseLevel1_flatStat_powerScale1.mat', 'dataFull_denoise')

%% Reconstruct the data

% in situ parameters
% Ni = 10000;
gamma = 0.5;
% plotFreq = inf;

static = objectTotal(:,:,1).*maskStatic;
dynamic = objectTotal.*maskDynamic;


[dynamicRec,staticRec,err] = inSituRecon_GPU_active(sqrt(dataFull_denoise),maskStatic,maskDynamic,static,probeTotal,[],Ni,probe.sphs0,gamma,plotFreq);

end
