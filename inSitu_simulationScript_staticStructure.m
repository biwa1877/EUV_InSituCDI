%% Simulation running script

%% Run simulations on power scaling for flat static region
tic;

noiseLevels = [0];%[0 1 2 3 4];
powerLevels = [0.5:0.5:10];
binning = 1;
structure = 'flat';
propDist = 0;
Ni = 5000;
plotFreq = inf;
whichAlgorithm = 'inSituUCLA';
whichDenoising = 'None';

filenameBase = ['JILAscript_edits_powerScan'];

Np = numel(powerLevels);
clear dynamicRecTotal staticRecTotal errTotal objectTotal maskDynamic errRealTotal maskStatic dynamicRecMin staticRecMin
for pp = 1:Np
    % Run the simulation and reconstruction
    [dynamicRecTotal{pp},staticRecTotal{pp},errTotal(:,pp),objectTotal{pp},maskDynamic{pp},maskStatic{pp}] = ...
        insitu_full_simulator(powerLevels(pp),structure,noiseLevels,propDist,binning,Ni,plotFreq,whichAlgorithm,whichDenoising);
    % Compare the final reconstruction to the object.
    errRealTotal(:,pp) = sqrt( mean( mean( abs( dynamicRecTotal{pp}(maskDynamic{pp}==1) - objectTotal{pp}(maskDynamic{pp}==1) ).^2 ) ) );
    
    % select the dynamic region of the reconstruction for saving
    [mDy,mDx] = find(maskDynamic{pp}==1);
    dRangeY = min(mDy(:)):max(mDy(:));
    dRangeX = min(mDx(:)):max(mDx(:));
    dynamicRecMin{pp} = dynamicRecTotal{pp}(dRangeY,dRangeX,:);
    
    % select the static region of the reconstruction for saving
    [mDy,mDx] = find(maskStatic{pp}==1);
    dRangeY = min(mDy(:)):max(mDy(:));
    dRangeX = min(mDx(:)):max(mDx(:));
    staticRecMin{pp} = staticRecTotal{pp}(dRangeY,dRangeX,:);
end



resultsFolder = 'simulationResults';
formatOut='_yyyymmdd_HHMMSS';
savefilename=strcat(filenameBase,'_',datestr(fix(clock),formatOut),'.mat');

things2Save = {'errTotal','structure','powerLevels','errRealTotal','dynamicRecMin','staticRecMin'};
save(fullfile(resultsFolder,savefilename),things2Save{:});

toc;

%% Plot the error

figure(1);
subplot(1,1,1);
bar(reordercats(categorical(structure),structure),errRealTotal)
ylabel('Final RMS error');

figure(2);
for pp = 1:Np
    semilogy(errTotal(:,pp),'lineWidth',3);
    hold on;
end
hold off;
legend(structure{:});
xlabel('Iteration');
ylabel('RMS error');