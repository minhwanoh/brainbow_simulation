function simulate_segmentImage

graphData.c                         = 2e-3;
graphData.colorRadiusForPure        = 20;
graphData.minSizeForPure            = 50;
graphData.maxPerim                  = 0.5;
graphData.spatialNhoodRadius        = sqrt(3)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.kMeansIterCount           = 1000;
disp(1)
cellCounts = [5 9 13];
rWSDs      = [0.01:0.01:0.04 0.1];
for incorrect = -3:3
  cC = 2;
  cellCount = cellCounts(cC);
  graphData.opts_irbleigs.K         = cellCount + incorrect;
  for channelCount = 3:5
    for rWSD = 1:numel(rWSDs)
      randomWalkSD = rWSDs(rWSD);
      for normalSD = 0.05:0.05:0.1
	 mergedSvFileName   = ['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2_sAff10.mat'];
         [index, graphData] = segmentImage(mergedSvFileName, graphData);
         load(mergedSvFileName, 'svCells');
         load(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_denoised.mat'], 'volumeLabels', 'denoised');
         clusterCount       = max(index);
         segmentation       = zeros(size(denoised,1),size(denoised,2),size(denoised,3)); for kka=1:clusterCount; thisCluster = find(index==kka); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = kka; end; end;
         allLabels=[]; allSeg=[];
         for kk = 1:numel(volumeLabels)
	   theseVox        = find(volumeLabels{kk});
           allLabels       = [allLabels; kk*ones(numel(theseVox),1)]; allSeg=[allSeg; segmentation(theseVox)];
         end
         toRemove          = find(allSeg==0);
         allLabels(toRemove)=[]; allSeg(toRemove)=[];
         [ARfg,RIfg,MI,HI] = valid_RandIndex(allLabels, allSeg);

         segmentationGT    = zeros(size(segmentation)); for kka=1:numel(volumeLabels); segmentationGT(volumeLabels{kka})=kka; end;
         segmentation(segmentation==0)     = max(segmentation(:))+1;
         segmentationGT(segmentationGT==0) = max(segmentationGT(:))+1;
         [AR,RI,MI,HI]     = valid_RandIndex(segmentation(:), segmentationGT(:));
save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_' num2str(graphData.opts_irbleigs.K) 'clusterCount_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
      end
    end
  end
end


%mergedSvFileName    = ['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
%thisCluster = find(index==0); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = clusterCount+1; end; % GARBAGE CLUSTER
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_merge1_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_garbage_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
