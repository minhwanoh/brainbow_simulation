clear all;

baseline_noises = [0, 200, 500];
protein_noises = [0, 0.1, 0.2];
protein_densities = [1e-6, 1e-7, 1e-8];
channelCounts = [3 4 5 6];
cellCounts = [5 9 13];

expansion_factor = '1.0';
baseline_noise = baseline_noises(1);
protein_noise = protein_noises(1);
protein_density = protein_densities(1);

for ch = 1:length(channelCounts)
    channelCount = channelCounts(ch);
    for cC = 1:length(cellCounts)
        cellCount = cellCounts(cC);
        load(['simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'])
        load(['simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2_sAff10.mat'])
        load(['simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat'])
        bbVol = overallRawVolume;
        bbVol(bbVol<0)=0;
        for kk = 1:size(bbVol, 4); 
            rawStack = bbVol(:,:,:,kk); 
            rawStack = rawStack - min(rawStack(:)); 
            rawStack = rawStack / max(rawStack(:)); 
            bbVol(:,:,:,kk) = rawStack; 
        end; 
        clear rawStack;
        clusterCount                                 = max(index);
        segmentation = zeros(stackSize); 
        for kka=1:clusterCount; 
            thisCluster = find(index==kka); 
            for kkb=1:numel(thisCluster); 
                segmentation(svCells{thisCluster(kkb)}) = kka; 
            end; 
        end;
        voxelCount                                   = numel(segmentation);
        xTileCount                                   = round(sqrt((clusterCount+1)/3) * 1) +2;
        yTileCount                                   = ceil((clusterCount+1)/xTileCount);
        compx                                        = size(segmentation,1)+size(segmentation,3);
        compy                                        = size(segmentation,2)+size(segmentation,3);
        bigIm                                        = ones( (size(segmentation,1)+1)*xTileCount-1, (size(segmentation,2)+1)*yTileCount-1, size(bbVol,4));
        xTile                                        = 1;
        yTile                                        = 1;
        bigIm(1:size(segmentation,1), 1:size(segmentation,2), :) = squeeze(max(bbVol, [], 3));
        for kk=1:clusterCount
          [xTile, yTile]                             = ind2sub([xTileCount yTileCount], kk+1);
          tmp2                                       = find(segmentation==kk);
          for mm = 1:3
            tmp1                                     = zeros(size(segmentation));
            tmp1(tmp2)                               = bbVol(tmp2+(mm-1)*voxelCount);
            bigIm((xTile-1)*(size(segmentation,1)+1)+1:xTile*(size(segmentation,1)+1)-1, (yTile-1)*(size(segmentation,2)+1)+1:yTile*(size(segmentation,2)+1)-1, mm) = max(tmp1, [], 3);
          end
        end
        imwrite(bigIm(:,:,1:3),['simData/bigim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef.png'])
    end
end