function simulate_segmentImage

graphData.c                         = 2e-3; % decay rate of color distancee
graphData.colorRadiusForPure        = 20; % threshold for max color dist
graphData.minSizeForPure            = 50; 
graphData.maxPerim                  = 0.5;
graphData.spatialNhoodRadius        = sqrt(3)+eps;
graphData.maxColorRadiusForProximal = 50;
graphData.minEdgeCountForProximal   = 5;
graphData.kMeansIterCount           = 1000;
%disp(1)

baseline_noises = [0, 200, 500];
protein_noises = [0, 0.1, 0.2];
protein_densities = [1e-6, 1e-7, 1e-8];
channelCounts = [3 4 5 6];
cellCounts = [5 9 13];

expansion_factor = '1.0';
baseline_noise = baseline_noises(1);
protein_noise = protein_noises(1);
protein_density = protein_densities(1);

AR_mat_cellCount = [];
AR_mat_baseline_noise = [];
AR_mat_protein_noise = [];
AR_mat_protein_density = [];


for ch = 1:length(channelCounts)
    channelCount = channelCounts(ch);
    baseline_noise = baseline_noises(1);
    protein_noise = protein_noises(1);
    protein_density = protein_densities(1);
    
    % cell count experiment
    for cC = 1:3
        cellCount = cellCounts(cC);
        [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, cellCount);
        ARfg_mat_cellCount(ch,cC) = ARfg; 
        AR_mat_cellCount(ch,cC) = AR;
    end

    % baseline_noise experiment
    cellCount = cellCounts(2);
    for bn = 1:length(baseline_noises)
        baseline_noise = baseline_noises(bn);
        [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, cellCount);
        ARfg_mat_baseline_noise(ch,bn) = ARfg; 
        AR_mat_baseline_noise(ch,bn) = AR;
    end
    
    % protein_noise experiment
    baseline_noise = baseline_noises(1);
    for pn = 1:length(protein_noises)
        protein_noise = protein_noises(pn);
        [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, cellCount);
        ARfg_mat_protein_noise(ch,pn) = ARfg; 
        AR_mat_protein_noise(ch,pn) = AR;
    end
    
    % protein density experiment
    protein_noise = protein_noises(1);
    for pd = 1:length(protein_densities)
        protein_density = protein_densities(pd);
        [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, cellCount);
        ARfg_mat_protein_density(ch,pd) = ARfg; 
        AR_mat_protein_density(ch,pd) = AR;
    end
    
        % protein density experiment
    protein_noise = protein_noises(1);
    for pd = 1:length(protein_densities)
        protein_density = protein_densities(pd);
        [ARfg, AR] = getSegmentResults(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor, graphData, cellCount);
        ARfg_mat_protein_density(ch,pd) = ARfg; 
        AR_mat_protein_density(ch,pd) = AR;
    end
    
end

t = datetime('now');
DateString = datestr(t);
save(['/vega/stats/users/mo2499/bbSimulation/simData/AR_results_' DateString '.mat'], 'ARfg_mat_cellCount', 'AR_mat_cellCount', 'ARfg_mat_baseline_noise', 'AR_mat_baseline_noise', 'ARfg_mat_protein_noise', 'AR_mat_protein_noise', 'ARfg_mat_protein_density', 'AR_mat_protein_density');
end

%mergedSvFileName    = ['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
%thisCluster = find(index==0); for kkb=1:numel(thisCluster); segmentation(svCells{thisCluster(kkb)}) = clusterCount+1; end; % GARBAGE CLUSTER
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_merge1_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
%save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_garbage_index.mat'], 'index', 'RIfg', 'ARfg', 'RI', 'AR', '-v7.3');
