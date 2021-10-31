function [dynamic,static,err] = inSituRecon_active(Intensity,maskStatic,maskDynamic,static,probe,dynamic,Ni)
%
% This function, INSITURECON_V3 runs the in-situ CDI algorithm created by Yuan
% Hung Lo and described in his ArXiv paper "In situ coherent diffractive 
% imaging" [https://arxiv.org/abs/1703.07219]
%
% This function uses a simplified algorithm that requires the probe and the
% static region to be known.
% 
% example usage:
% [Intensity,maskStatic,maskDynamic,dynamic,static,probe] = inSituSimulation_v3();
% [dynamicRecovered,staticRecovered,error] = inSituRecon_v3(Intensity,maskStatic,maskDynamic,static,probe,[]);
% min(error(:))<1e-10
%
% The above three lines should return 1. If they return 0, something has
% gone awry.
%
%
% inputs:
%
%       Intensity                   An Ny x Nx x Nt array containing the
%                                   diffraction patterns in intensity (not
%                                   amplitude)
%
%       maskStatic                  An Ny x Nx array identifying the
%                                   regions that do not change frame to
%                                   frame.
%
%       maskDynamic                 An Ny x Nx array identifying the
%                                   regions that can change from frame to
%                                   frame. Note: maskDynamic+maskStatic
%                                   should be the entire FOV on the sample.
%
%       static                      Initital guess for the static region of
%                                   the object.
%
%       probe                       Initial guess for the probe.
%
%       dynamic                     Initial guess for the dynamic region of
%                                   the object.
%
%       Ni                          Number of iterations to run.
%
%
% outputs:
%
%       dynamic                     Recontructed dynamic region.
%
%       static                      Reconstructed static region.
%
%       err                         Error as a function of iteration.
%
%
% Version 1: first stable version. Note that after 250 iterations, this
% version should return machine precision results. 
%
% Version 2: includes a probe. If the correct probe is given, the algorithm
% reaches machine precsion after ~3700 iterations.
%
% Version 3: uses a larger gaussian, complex probe and a complex, random 
% object. Reaches machine precsion after 250 iterations. Added option to
% feed in a guess for the dynamic portion.
%
% See also:
%
% INSITUSIMULATION_V3
%
% Author information:
%
%       Robert Karl: robert.karl@colorado.edu
%
% 03-08-2018
%

%% Settings

% Ni = 10; %number of iterations
plotFreq = 5;

%% input parsing

[Ny,Nx,Nt] = size(Intensity);
% Take the square root!
ESWft = sqrt(Intensity);

%% Create initial guess

static = repmat(static,[1,1,Nt]);
if isempty(dynamic)
    dynamic = repmat(maskDynamic,[1,1,Nt]);
end

%% The phase retrieval part

gamma = 0.5; %time relaxation factor

%Initialize figures for plotting
figure(1);close(1);figure(1);
figure(2);close(2);figure(2);
figure(3);close(3);figure(3);

% Determine dynamic region to only plot there
[mDy,mDx] = find(maskDynamic==1);
dRangeY = min(mDy(:)):max(mDy(:));
dRangeX = min(mDx(:)):max(mDx(:));

% preallocate
StP = zeros(Ny,Nx,Nt);
object = zeros(Ny,Nx,Nt);
psi = zeros(Ny,Nx,Nt);
Psi = zeros(Ny,Nx,Nt);
PsiP = zeros(Ny,Nx,Nt);
psiP = zeros(Ny,Nx,Nt);
objectP = zeros(Ny,Nx,Nt);

err = zeros(Ni,1);

eps = 1e-9;

        % Calculating the probe ratio
        Pratio = conj(probe)./(abs(probe)+eps).^2;

for ii = 1:Ni
    for tt = 1:Nt
        temp = circshift(static,[0,0,1]);
        StP(:,:,tt) = gamma*temp(:,:,tt)+(1-gamma)*static(:,:,tt);
        object(:,:,tt) = StP(:,:,tt)+dynamic(:,:,tt);
        % multiply by probe
        psi(:,:,tt) = object(:,:,tt).*probe;
        Psi(:,:,tt) = fftshift(fft2(ifftshift(psi(:,:,tt))));
        % Apply Fourier constraint
        PsiP(:,:,tt) = ESWft(:,:,tt).*Psi(:,:,tt)./abs(Psi(:,:,tt));
        psiP(:,:,tt) = fftshift(ifft2(ifftshift(PsiP(:,:,tt))));
        psiDel = psiP(:,:,tt)-psi(:,:,tt);


        objectP(:,:,tt) = object(:,:,tt) + Pratio.*psiDel;
        
        static(:,:,tt) = objectP(:,:,tt).*maskStatic;
        dynamic(:,:,tt) = objectP(:,:,tt).*maskDynamic;
    
    end
    
    % Calculate error: normalized by # of pixels
    err(ii) = real(sqrt(mean( (abs(Psi(:))-ESWft(:)).^2 ))/(Ny*Nx));
    
    % Plot the status
    if mod(ii,plotFreq)==0
        % Plot the magnitude
%         figure(1);
%         for tt = 1:Nt
%             subplot(ceil(sqrt(Nt)),ceil(sqrt(Nt)),tt);
%             imagesc(abs(dynamic(dRangeY,dRangeX,tt)));
%             title(['frame #',num2str(tt)]);
%         end
%         drawnow;
%         % Plot the phase
%         figure(2);
%         for tt = 1:Nt
%             subplot(ceil(sqrt(Nt)),ceil(sqrt(Nt)),tt);
%             imagesc(angle(dynamic(dRangeY,dRangeX,tt)));
%             title(['frame #',num2str(tt)]);
%         end
%         drawnow;
        % Plot the error
        figure(3);
        semilogy(err(1:ii));
        drawnow;
    end

    
end

