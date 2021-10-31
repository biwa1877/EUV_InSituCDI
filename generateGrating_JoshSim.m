%% Simulating Acoustic dynamics for in-situ CDI
%
% This function takes Josh Knobloch's COMSOL simulations on a nickel line
% and prepares them for in-situ simulations. 
% The output of this function is the object and probe for the in-situ
% experiment.
%
% Outputs:
%
%           objectTotal             Npix x Npix x Nf array with the full
%                                   object (static and dynamic) for each
%                                   time frace (Nf total).
%
%           probeTotal              Npix x Npix array with the full probe
%                                   (static and dynamic).
%
%           maskStatic              Npix x Npix array where the static
%                                   region is 1 and the rest is 0.
%
%           maskDynamic             Npix x Npix array where the dynamic
%                                   region is 1 and the rest is 0.
%

function [objectTotal,probeTotal,maskStatic,maskDynamic] = generateGrating_JoshSim()

%% Simulation vs. time Trace

% Load the simulation:
longPath = 'JoshSimulation';
% 2.2um wide line
% filepath = 'DynamicImaging_FlagSimulation_JelloModelT450_2DBiggerSlice_FinerMesh3_0by0pt1to600ps_Good_Surf.txt';
% 1.5 um wide line
filepath = 'FlagSimulation_JelloModelT450K_Mesh3_W1pt5um_H20pt5nm_T0by0pt1to600ps_Surface.txt';
dt = 0.1*10^-12; %delta time for the simulations
dy = 6e-9; % delta space for the simulations

%now let us read it into a data array.
data = importdata(fullfile(longPath,filepath));

%get some dimensions
[Np,Num] = size(data);
Nt = (Num-2)/2; %the number of time positions

%now sort it
xposc = data(:,1); %x position
zposc = data(:,2); %axial position

xdispc = zeros(Np,Nt); %x displacement
zdispc = zeros(Np,Nt); %axial displacement
for i = 1:Nt
    xdispc(:,i) = data(:,2*i+1);
    zdispc(:,i) = data(:,2*i+2);
end

%now build the time matrix
time = (0:1:Nt-1)*dt; %in ps for

% dealing with troublesome points
% These two points on the edge are out of order, leading to an awkward
% discontinuity, but if I just switch their order all should be well
swap1 = 352;%469;
swap2 = 353;%470;

temp = xdispc(swap1,:);
xdispc(swap1,:) = xdispc(swap2,:);
xdispc(swap2,:) = temp;
temp = zdispc(swap1,:);
zdispc(swap1,:) = zdispc(swap2,:);
zdispc(swap2,:) = temp;
temp = xposc(swap1,:);
xposc(swap1,:) = xposc(swap2,:);
xposc(swap2,:) = temp;
temp = zposc(swap1,:);
zposc(swap1,:) = zposc(swap2,:);
zposc(swap2,:) = temp;

%% isolate regions that are nickel and silicon

maskNickel = zeros(size(zposc));

nickelStart = 102;
nickelEnd = 352;
maskNickel(nickelStart:nickelEnd,:) = 1;
maskSilicon = 1-maskNickel;

%% Select time steps to use

Nf = 9;
timeSpacing = 50E-12;

if (Nf-1)*(timeSpacing/dt)+1 > Nt
    error("You've asked to go beyond the simulation window. Either use fewer frames (Nf) or smaller time spacing");
end

timesUse = [0:Nf-1]*(timeSpacing/dt);
timesUse = round(timesUse)+1;

zdispc_frames = zdispc(:,timesUse);
time_frames = time(timesUse);

%% Form the full object

lambda = 28.9E-9;
theta = 45*pi/180;

% Calculate geometric phase
height = maskNickel*20e-9 + zdispc_frames;
geometricPhase = 4*pi*cos(theta)*height/lambda;

n_Ni = 1-0.1386720240 - 1i*0.111483328; % From CXRO
n_Si = 1-0.0628686249 - 1i*0.008287034; % From CXRO

% Calculate Fresnel reflectivities
cosThNi = sqrt(1-n_Ni^2*sin(theta)^2);
cosThSi = sqrt(1-n_Si^2*sin(theta)^2);

r_Ni = ( n_Ni*cos(theta)-cosThNi )./( n_Ni*cos(theta)+cosThNi );
r_Si = ( n_Si*cos(theta)-cosThSi )./( n_Si*cos(theta)+cosThSi );

% Combine for full complex reflectivity
object1D = r_Ni.*maskNickel + r_Si.*maskSilicon;
% object1D = object1D.*(1+zdispc_frames/max(zdispc_frames(:))*0.3);
object1D = object1D.*exp(1i.*geometricPhase);

%% Create the rest of the object so that we have a static region

% Experimental parameters
D2Det = 4E-2;
detPix = 13.5e-6;
Npix = 2048;

pix = lambda*D2Det/(Npix*detPix);

%% Extrude the 1D simulation to make a series of 2D images

[Ny,~] = size(object1D);

object2D = zeros(Ny,Ny,Nf);
for ff = 1:Nf
    object2D(:,:,ff) = repmat(object1D(:,ff),[1,Ny]);
end

%% Interpolating onto the experimental reconstruction grid
[Ny,Nx,Nf] = size(object2D);
gridSimY = (1:Ny)*dy;
gridSimX = (1:Nx)*dy;
[gridSimX,gridSimY] = meshgrid(gridSimX,gridSimY);

gridExpY = pix:pix:Ny*dy;
gridExpX = pix:pix:Nx*dy;
[gridExpX,gridExpY] = meshgrid(gridExpX,gridExpY);

[Nyo,Nxo] = size(gridExpX);
object2DExp = zeros(Nyo,Nxo,Nf);
for ff = 1:Nf
    object2DExp(:,:,ff) = interp2(gridSimX,gridSimY,object2D(:,:,ff),gridExpX,gridExpY);
end

% Repeat the simulation to make a full grating
NrepY = 5;
NrepX = 30;
object2DExpRep = zeros(NrepY*Nyo,NrepX*Nxo,Nf);
for ff = 1:Nf
    object2DExpRep(:,:,ff) = repmat(object2DExp(:,:,ff),[NrepY,NrepX]);
end

% Make the rest of the object be substrate (Only effective if interpolation
% grids are mismatched)
object2DExpRep(isnan(object2DExpRep)) = r_Si;

[Ny,Nx,Nf] = size(object2DExpRep);
% pad the array so that we have the same number of pixels as the detector.
object2DExpRep = padarray(object2DExpRep,[Npix-Ny,Npix-Nx],r_Si,'post');

% Shift the object so that the dynamic region is not on the edge
for ff = 1:Nf
    object2DExpRep(:,:,ff) = circshift(object2DExpRep(:,:,ff),floor([Npix/2,0]));
end

%% Create the probes
% Using a split mirror geometry, we will have a hard edge on one side of
% the probe, on each probe. This will be the interior edge of the probe.


beamWidth = 40e-6;
hardEdgeWidth = 7E-6;
sigma = lambda*D2Det/(detPix*hardEdgeWidth);

[pX,pY] = meshgrid((-Npix/2:Npix/2-1)*pix,(-Npix/2:Npix/2-1)*pix);

probe = exp(-(pX/beamWidth).^2-(pY/beamWidth).^2);
% square support
% probe(abs(pX)>beamWidth | abs(pY)>beamWidth)=0;
% circle support
probe( (abs(pX).^2+abs(pY).^2) > hardEdgeWidth.^2) = 0;
% Cut the beam in half
probeStatic = probe;
probeStatic(pY>0) = 0;
probeStatic = circshift(probeStatic,[-floor(hardEdgeWidth/2/pix),0]);

probeDynamic = probe;
probeDynamic(pY<0) = 0;

probeTotal = probeStatic+probeDynamic;

%% Create the static and dynamic masks

maskStatic = probeStatic>0;
maskDynamic = probeDynamic>0;

%% Parse outputs
% This section renames variables so that they are user-friendly outside of
% this script. It also serves as a handy reference for what was generated.

objectTotal = object2DExpRep;
probeTotal = probeTotal;
maskStatic = maskStatic;
maskDynamic = maskDynamic;

end