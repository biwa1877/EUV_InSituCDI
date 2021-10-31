%% Simulation running script

%% Run simulations on different noise level
tic;

noiseLevels = [1];%[0 1 2 3 4];
Nn = numel(noiseLevels);

powerLevels = [0.01 0.1 1 10];   %[0.001 0.01 0.1 1 10];
Np = numel(powerLevels);

binning = 1;

% static structure:
%   could be:   'flat'
%               'speckle'
%               'cameraman'
%               'vertGrating'
%               'horGrating'
structure = 'speckle';

propDist = 0;   % probe propagation distance from the pinhole, 1e-6;
Ni = 10;     % max #iteration
plotFreq = 10;  % plot frequency

% 'inSitu' or 'inSituRAAR'
whichAlgorithm = 'inSituRAAR';

% 'BM4D' or 'Thresholding' or 'None'
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%           BM4D on 2048x2048x9 dataset will take ~24min
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
whichDenoising = 'BM4D';

for pp = 1:Np   % power scan loop
    for nn = 1:Nn   % noise scan loop
        % house keeping
        clear dynamicRecTotal staticRecTotal errTotal objectTotal maskDynamic errRealTotal maskStatic dynamicRecMin staticRecMin
        
        % Run the simulation and reconstruction
        [dynamicRecTotal{nn},staticRecTotal{nn},errTotal(:,nn),objectTotal{nn},maskDynamic{nn},maskStatic{nn}] ...
            = insitu_full_simulator(powerLevels(pp),structure,noiseLevels(nn),propDist,binning,Ni,plotFreq,whichAlgorithm,whichDenoising);
        % Compare the final reconstruction to the object.
        errRealTotal(:,nn) = sqrt( mean( mean( abs( dynamicRecTotal{nn}(maskDynamic{nn}==1) - objectTotal{nn}(maskDynamic{nn}==1) ).^2 ) ) );
        
        % select the dynamic region of the reconstruction for saving
        [mDy,mDx] = find(maskDynamic{nn}==1);
        dRangeY = min(mDy(:)):max(mDy(:));
        dRangeX = min(mDx(:)):max(mDx(:));
        dynamicRecMin{nn} = dynamicRecTotal{nn}(dRangeY,dRangeX,:);
        
        % select the static region of the reconstruction for saving
        [mDy,mDx] = find(maskStatic{nn}==1);
        dRangeY = min(mDy(:)):max(mDy(:));
        dRangeX = min(mDx(:)):max(mDx(:));
        staticRecMin{nn} = staticRecTotal{nn}(dRangeY,dRangeX,:);
    end
    
    resultsFolder = 'InSituCDI with Noise/Speckle As Static Structure/BM4D Denoising';
    % filenameBase = 'noise_flatStatic';
    filenameBase = ['noisy_',structure,'StaticStruct_',num2str(powerLevels(pp)),'PowerScale'];
    formatOut='yyyymmdd_HHMMSS';
    savefilename=strcat(filenameBase,'_',datestr(fix(clock),formatOut),'.mat');
    
    things2Save = {'errTotal','structure','noiseLevels','errRealTotal','dynamicRecMin','staticRecMin'};
    save(fullfile(resultsFolder,savefilename),things2Save{:});
end
toc;

%% Plot the error

figure(1);
plot(powerLevels,errRealTotal);

figure(2);
surf(powerLevels,1:Ni,log(errTotal),'edgeColor','interp');
