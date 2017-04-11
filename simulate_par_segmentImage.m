function simulate_par_segmentImage(paramID)

graphData.c                         = 2e-3; % decay rate of color distancee
graphData.colorRadiusForPure        = 20; % threshold for max color dist
graphData.minSizeForPure            = 50; 
graphData.maxPerim                  = 0.5;
graphData.spatialNhoodRadius        = sqrt(3)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.kMeansIterCount           = 1000;
%disp(1)

channelCounts = [3 4 5 6];
cellCounts = [5 9 13];
baseline_noises = [0, 200, 500];
protein_noises = [0, 0.1, 0.2];
protein_densities = [1e-6, 1e-7, 1e-8];

ind = mod(paramID, 12);
if (ind == 0)
    ind = 12;  paramID = paramID - 12;
end
ii = floor(paramID/12);
[kk1 kk2] =ind2sub([numel(channelCounts), 3], ind);
paramIndx = [1 2 1 1 1];
paramIndx(1) = kk1; 
paramIndx(ii+2) = kk2;

expansion_factor = '1.0';
channelCount = channelCounts(paramIndx(1));
cellCount = cellCounts(paramIndx(2));
baseline_noise = baseline_noises(paramIndx(3));
protein_noise = protein_noises(paramIndx(4));
protein_density = protein_densities(paramIndx(5));

load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'volumeLabels', 'denoised');
actualCellCount = length(find(~cellfun(@isempty,volumeLabels)))
graphData.opts_irbleigs.K         = actualCellCount; 
mergedSvFileName   = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2_sAff10.mat'];
[index, graphData] = segmentImage(mergedSvFileName, graphData);
load(mergedSvFileName, 'svCells');

fground = 0;
for i = 1:length(volumeLabels)
    fground = fground + nnz(volumeLabels{i}(:,:,25:75)); 
end
stackSize = size(denoised);
expressionDensity = fground/(stackSize(1)*stackSize(2)*50);

intersectingVox=[]; 
for kk1=1:numel(volumeLabels)
    for kk2=kk1+1:numel(volumeLabels)
        intersectingVox = [intersectingVox; intersect(find(volumeLabels{kk1}),find(volumeLabels{kk2}))];
    end;
end;

seg1=zeros(size(denoised,1),size(denoised,2),size(denoised,3));
seg2=zeros(size(denoised,1),size(denoised,2),size(denoised,3));
for kk=1:numel(svCells)
    seg1(svCells{kk})=index(kk);
end;
for kk=1:numel(volumeLabels)
    seg2(volumeLabels{kk})=kk;
end;
active=intersect(find(seg1),find(seg2));
active=setdiff(active,intersectingVox);

tt1=seg1(active); 
tt2=seg2(active);
[ARfg,RIfg,MI,HI] = valid_RandIndex(tt1,tt2);

segmentationGT    = zeros(size(seg1)); 
for kka=1:numel(volumeLabels); 
    segmentationGT(volumeLabels{kka})=kka;
end;
seg1(seg1==0)     = max(seg1(:))+1;
segmentationGT(segmentationGT==0) = max(segmentationGT(:))+1;
[AR,RI,MI,HI]     = valid_RandIndex(seg1(:), segmentationGT(:));
save(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3', 'expressionDensity', 'actualCellCount');

end