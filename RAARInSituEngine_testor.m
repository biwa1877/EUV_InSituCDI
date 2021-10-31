%% insitu Noise simulator
% This function take the obejct and probe from the generator function and
% creates diffraction patterns with realistic noise!! (or no noise if you
% want that).

function [dynamicRec,staticRec,err] = RAARInSituEngine_testor(structure,powerScale)

%% default inputs
if ~exist('structure','var')
    structure = 'flat';
end
if ~exist('powerScaling','var')
    powerScale = 10;
end

%% Add required paths
addpath('Ptychography Data Generation Code');

%% Generate the object and probe
[objectTotal,probeTotal,maskStatic,maskDynamic] = generateGrating_JoshSim();

%% Change the static structure

% structure = 'flat'; % Options are 'flat' 'speckle' and 'cameraman'

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
end

probeTotal(maskStatic==1) = probeTotal(maskStatic==1).*temp(maskStatic==1);

%% Change the power scaling

powerStatic = sum(sum(abs(probeTotal.*objectTotal(:,:,1).*maskStatic).^2));
[~,~,Nf] = size(objectTotal);

powerDynamic = zeros(1,Nf);
for ff = 1:Nf
    powerDynamic(ff) = sum(sum(abs(probeTotal.*objectTotal(:,:,ff).*maskDynamic).^2));
end

probeTotal(maskStatic==1) = probeTotal(maskStatic==1).*mean(powerDynamic)/powerStatic;


% powerScale = 10;
probeTotal(maskStatic==1) = powerScale.*probeTotal(maskStatic==1);

%% Propagate probe slightly


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

camSetting.pnsFlag = 0; %Poisson noise
camSetting.wgnFlag = 0; % White gaussian noise
camSetting.readRate = 0; % Readout rate noise: 0-4 depending on which rate
camSetting.N = Ny;
camSetting.Dt_exp = 1; %seconds?
camSetting.numAcc = 3;

%% Generate a diffraction pattern for each time frame

[Ny,Nx,Nf] = size(objectTotal);
dataFull = zeros(Ny,Nx,Nf);

for ff = 1:Nf

    % =================== begin object ===================

    object.obj0 = objectTotal(:,:,ff);

    % generate diffraction data
    [camSetting, probe, object, scanGrid, data] = ...
        GenPtyScanFFTShift_v01(camSetting, probe,object, scanGrid, data);
    
    dataFull(:,:,ff) = data.I;
end

%% change the scaling of the data
% vacuum permittivity
eps0 = 8.8542e-12; % A^2 * s^4/kg/m^3 OR F/m

% vacuum permeability
mu0 = 1.2566e-6; % N/A^2

% speed of light
c = sqrt(1/eps0/mu0); % m/s


% undoing the physical scaing from the realistic simulation to match
% propagator in the reconstruction code.

% Factors from the probe normalization (is squared later because these were
% imposed in amplitude)
probeNorm = sqrt( 2*camSetting.hw/(c*eps0) * probe.nPh_s/sum( abs(probeTotal(:)).^2 * probe.dxs * probe.dys ) );
% Factors from the ESW propagator
eswNorm = (c*eps0/2)/probe.lambda^2/probe.D2D^2*probe.dxs^2*probe.dys^2*(camSetting.N^2);
% Factors from the camera exposure
camNorm = camSetting.QE*camSetting.Dt_exp*camSetting.dxd^2/camSetting.hw;

dataFull = dataFull/(probeNorm^2)/eswNorm/camNorm; 
% dataFull = dataFull/eswNorm;

% %% diffraction patterns skipping noise simulation
% 
% norm = sqrt(camSetting.N^2);
% for ff = 1:Nf
%     dataFull(:,:,ff) = objectTotal(:,:,ff).*probeTotal.*probe.sphs0;
%     dataFull(:,:,ff) = fftshift(fft2(ifftshift(dataFull(:,:,ff))))/norm;
%     dataFull(:,:,ff) = abs(dataFull(:,:,ff)).^2;
% end

%% Reconstruct the data

% in situ parameters
Ni = 5000;
gamma = 0.5;
plotFreq = 10;  %inf;

% RAAR parameters
beta = 1;    % step size
eta1 = 0.8;    % modulus reflection depth
eta2 = -1;    % overlap reflection depth

static = objectTotal.*maskStatic;
dynamic = objectTotal.*maskDynamic;

% [dynamicRec,staticRec,err] = inSituRecon_GPU_active(sqrt(dataFull),maskStatic,maskDynamic,static,probeTotal,[],Ni,probe.sphs0,gamma,plotFreq);
[dynamicRec,staticRec,err] = inSitu_RAAR_Recon_GPU_active_v04(sqrt(dataFull),...
    maskStatic,maskDynamic,static,probeTotal,[],Ni,probe.sphs0,gamma,beta,eta1,eta2,plotFreq);

end
