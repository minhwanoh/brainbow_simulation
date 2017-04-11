function simulate_segmentImage_1

graphData.c                         = 2e-3; % decay rate of color distancee
graphData.colorRadiusForPure        = 20; % threshold for max color dist
graphData.minSizeForPure            = 50; 
graphData.maxPerim                  = 0.5;
graphData.spatialNhoodRadius        = sqrt(3)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.kMeansIterCount           = 1000;
disp(1)

cellCounts = [5 9 13]; % [2 5 9 13];
actualCellCounts = [4 9 11];
for cC = 1:3 % 4
    cellCount = cellCounts(cC);
   % ActualcellCount = actualCellCounts(cC);
    load(['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sim_' num2str(cellCount) 'cells_ch4_denoised.mat'], 'volumeLabels', 'denoised');
    graphData.opts_irbleigs.K         = length(find(~cellfun(@isempty,volumeLabels))); % 
    mergedSvFileName   = ['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sim_' num2str(cellCount) 'cells_ch4_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2_sAff10.mat'];
    [index, graphData] = segmentImage(mergedSvFileName, graphData);
    load(mergedSvFileName, 'svCells');

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
    save(['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sim_' num2str(cellCount) '_ch4_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
end 
end

%mergedSvFileName    = ['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
%thisCluster = find(index==0); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = clusterCount+1; end; % GARBAGE CLUSTER
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_merge1_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_garbage_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
