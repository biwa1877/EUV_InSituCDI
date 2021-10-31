%% this version 03 is the frist version that works and works better than the original insitu code;
%   But there is still a problem with it, which is that it only works for
%   eta = -1. There has something to do with the cross talk between the
%   RORM_psi(t) and psi(t+1). I will try to address this issue in version
%   04.

function [dynamic,static,err] = inSitu_RAAR_Recon_GPU_active_v03(ESWft,maskStatic,maskDynamic,static,probe,dynamic,Ni,sphs0,gamma,beta,eta1,eta2,plotFreq)
%
% This function, INSITURECON_V3 runs the in-situ CDI algorithm created by Yuan
% Hung Lo and described in his ArXiv paper "In situ coherent diffractive 
% imaging" [https://arxiv.org/abs/1703.07219]
%
% This function uses a simplified algorithm that requires the probe and the
% static region to be known.
%
% This version runs on a GPU.
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
%       ESWft                       An Ny x Nx x Nt array containing the
%                                   diffraction patterns in amplitude (not
%                                   intensity)
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
%       sphs0                       Fresnel phase for propagation.
%
%       gamma                       Relaxation factor for static region
%                                   between times
%
%       plotFreq                    Number of iterations before the plot
%                                   updates. Set to inf to turn off
%                                   plotting.
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
% plotFreq = inf;

%% input parsing

[Ny,Nx,Nt] = size(ESWft);
ESWft = fftshift(fftshift(ESWft,1),2);
norm = sqrt(Ny*Nx);
% Take the square root!
% ESWft = sqrt(Intensity);
maskDynamic = double(maskDynamic);
maskStatic = double(maskStatic);

%% Create initial guess

static = repmat(double(static),[1,1,Nt]);
if isempty(dynamic)
    dynamic = repmat(maskDynamic,[1,1,Nt]);
end

%% The phase retrieval part

% gamma = 0.8; %time relaxation factor

%Initialize figures for plotting
% figure(1);close(1);figure(1);
% figure(2);close(2);figure(2);
% figure(3);close(3);figure(3);

% Determine dynamic region to only plot there
[mDy,mDx] = find(maskDynamic==1);
dRangeY = min(mDy(:)):max(mDy(:));
dRangeX = min(mDx(:)):max(mDx(:));

% preallocate
% StP = zeros(Ny,Nx,Nt);
% object = zeros(Ny,Nx,Nt);
% psi = zeros(Ny,Nx,Nt);
% Psi = zeros(Ny,Nx,Nt);
% PsiP = zeros(Ny,Nx,Nt);
% psiP = zeros(Ny,Nx,Nt);
% objectP = zeros(Ny,Nx,Nt);

err = zeros(Ni,1);

eps = 1e-14;

% Calculating the probe ratio
Pratio = abs(probe) / max(abs(probe(:))) .* conj(probe)./(abs(probe)+eps).^2;
        
% Put things onto the GPU
[static,gamma,dynamic,probe,ESWft,Pratio,maskStatic,maskDynamic,Ny,Nx,err,sphs0] = ...
    send2GPU(static,gamma,dynamic,probe,ESWft,Pratio,maskStatic,maskDynamic,Ny,Nx,err,sphs0);

% Construct the figure
        figure(1);close(1);figure(1);
        subplot(121);
%         imagesc(abs(dynamic(dRangeY,dRangeX,tt).*probe(dRangeY,dRangeX)));
        iAx = imagesc(abs(dynamic(dRangeY,dRangeX,1)));
        subplot(122);
        eAx = semilogy(err(1));
        drawnow;

% RAAR parameter
if isempty(beta)
    beta = gpuArray(0.75);
else
    beta = gpuArray(beta);
end
if isempty(eta1)
    eta1 = gpuArray(eta1);
else
    eta1 = gpuArray(eta1);
end
if isempty(eta2)
    eta2 = gpuArray(0.75);
else
    eta2 = gpuArray(eta2);
end

for ii = 1:Ni
    for tt = 1:Nt
%         temp = circshift(static,[0,0,1]);
        StP = gamma * static(:,:,mod(tt+(Nt-2),Nt)+1)...
            + (1-gamma) * static(:,:,tt);
            % the expression "mod(tt+(Nt-2),Nt)+1" produces tt-1 for every
            % value, but for tt=1 it produces Nt instead.
%         StP = static(:,:,1);
        object = StP + dynamic(:,:,tt);
        
        % make ESW
        psi = object .* probe;
        
        if ii > 1
            % overlap reflector
            psi_OR = (1 + eta2) * psi_OP - eta2 * psi;
            % RAAR
            psi = 0.5 * beta * (psi_OR + psi) + (1 - beta) * psiP;
        end
        
%         Psi = fftshift(fft2(ifftshift(sphs0.*psi)))/norm;
        Psi = fft2(sphs0 .* psi) / norm;
        data = ESWft(:,:,tt);
        % Calculate error: normalized by # of pixels
        err(ii) = err(ii)+real(sqrt(mean( (abs(norm2max(Psi(:)))-norm2max(data(:))).^2 ))/(Ny*Nx));
        % Apply Fourier constraint
        PsiP = Psi;
        PsiP(data~=-1) = data(data~=-1).*exp(1i*angle(Psi(data~=-1)));
%         psiP = fftshift(ifft2(ifftshift(PsiP))).*conj(sphs0)*norm;
        psiP = ifft2(PsiP) .* conj(sphs0) * norm;
        
        psi_RM = (1 + eta1) * psiP - eta1 * psi;
%         psiDel = psiP-psi;
        psiDel = psi_RM - psi;
        
        objectP = object + Pratio .* psiDel;
        
        static(:,:,tt) = objectP .* maskStatic;
        dynamic(:,:,tt) = objectP .* maskDynamic;
        
        % overlap projector
        psi_OP = (static(:,:,tt) + dynamic(:,:,tt)) .* probe;
        
    end
    

    
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
%         figure(3);
        subplot(121);
%         imagesc(abs(dynamic(dRangeY,dRangeX,tt).*probe(dRangeY,dRangeX)));
%         imagesc(abs(dynamic(dRangeY,dRangeX,tt)));
        set(iAx,'CData',gather(abs(dynamic(dRangeY,dRangeX,mod(ii-1,Nt)+1))));
        title(num2str(mod(ii-1,Nt)+1));
        subplot(122);
        semilogy(err(1:ii));axis tight;
        drawnow;
    end

    
end

% Gather the outputs from the GPU
[dynamic,static,err] = getFromGPU(dynamic,static,err);


% End the main function
end


function [varargout] = send2GPU(varargin)
    for aa = 1:nargin
        varargout{aa} = gpuArray(varargin{aa});
    end
end

function [varargout] = getFromGPU(varargin)
    for aa = 1:nargin
        varargout{aa} = gather(varargin{aa});
    end
end

function a = norm2max(b)
% normalize input to maximum
    a = ( b - min(abs(b(:))) )/( max(abs(b(:))) - min(abs(b(:))) );
end