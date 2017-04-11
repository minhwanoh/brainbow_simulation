function [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, clusterSize)

load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'volumeLabels', 'denoised');
actualCellCount = length(find(~cellfun(@isempty,volumeLabels)))
graphData.opts_irbleigs.K         = clusterSize;
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
