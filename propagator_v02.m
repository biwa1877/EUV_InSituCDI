%% wave propagator

% add kz as an input option for birefringent crystals


function E2 = propagator_v02(E1,dx,lambda,z,n,kz)
    
    % default refractive index n
    if isempty(n) == 1
        n = 1;
    end
    
    % number of sampling points
    N = size(E1,1);
    
    % create a structure of parameters necessary for Parseval's theorem
    paras.dx = dx;
    paras.lambda = lambda;
    paras.z = z;
    
    % FFT
    [F1, dfx] = nfft_v02(E1, paras);
    
    % transfer function
    if isempty(kz)
        fx = (-N/2:N/2-1)' / (N * dx);
        [fx, fy] = meshgrid(fx);
        H = exp(1i * abs(2 * pi * sqrt((n/lambda).^2-fx.^2-fy.^2) * z));
    else
        H = exp(1i * kz * z);
    end
    % transfer function
    F2 = F1 .* H;
    
    % create a structure of parameters necessary for Parseval's theorem
    paras.dx = dfx;
    
    % iFFT
    E2 = nifft_v02(F2, paras);
    
end

function [F1, dfx] = nfft_v02(E1, paras)
    F1 = ifftshift(fftn(fftshift(E1)));
    dfx = 1;
end

function E2 = nifft_v02(F2, paras)
    E2 = ifftshift(ifftn(fftshift(F2)));
end