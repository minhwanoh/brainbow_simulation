function simulate_mergeSupervoxels_subgraph

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
    chCount = channelCount;
    mergeOpts.zAnisotropy                                                    = 5/4;
    mergeOpts.demix.maxSimilarNeighborNormLUVDist                            = 50 * sqrt(chCount/4);
    mergeOpts.demix.minImprovementFactor                                     = 5;
    mergeOpts.demix.maxSizeForDemixing                                       = 500;
    mergeOpts.mergeSmallSuperVoxels.luvColorDistanceUpperBound               = 20;
    mergeOpts.mergeSmallSuperVoxels.disconnectedSVsizeTh                     = 20;
    mergeOpts.mergeSmallSuperVoxels.maxVoxColorDist                          = 0.25; % 0.5;
    mergeOpts.mergeWRTnAo.sDist                                              = sqrt(3);
    mergeOpts.mergeWRTnAo.minDotProduct                                      = 0.9659; % pi/12 % sqrt(3)/2;
    mergeOpts.mergeWRTnAo.maxColorDist                                       = 10 * sqrt(chCount/4);
    mergeOpts.mergeWRTnAo.normFlag                                           = true;
    mergeOpts.mergeWRTnAo.maxVoxColorDist                                    = 0.25; % 0.5;
    mergeOpts.mergeWRTnAo.zAnisotropy                                        = mergeOpts.zAnisotropy;
    mergeOpts.mergeSingleNeighborSuperVoxels.maxVoxColorDist                 = 0.25; % 0.5;
    mergeOpts.mergeSingleNeighborSuperVoxels.maxSizeForSingleNeighborSVs     = 60;
    mergeOpts.mergeCloseNeighborhoods.maxDistNormLUV                         = 10 * sqrt(chCount/4);
    mergeOpts.mergeCloseNeighborhoods.maxVoxColorDist                        = 0.25; % 0.5;
    mergeOpts.kmeansMerging.overSamplingFactor                               = 5;
    mergeOpts.kmeansMerging.maxColorDistance                                 = 15 * sqrt(chCount/4);
    mergeOpts.kmeansMerging.maxSubgraphSize                                  = 10000;
    mergeOpts.kmeansMerging.cutSubgraphsMaxLuvColorDistance                  = 60;
    mergeOpts.spatialDistanceCalculation.upperBound                          = 2;
    mergeOpts2                                                               = mergeOpts;
    mergeOpts2.kmeansMerging.overSamplingFactor                              = 5;
    mergeOpts2.spatialDistanceCalculation.upperBound                         = 10;

    baseline_noise = baseline_noises(1);
    protein_noise = protein_noises(1);
    protein_density = protein_densities(1);
    
    for cC = 1:3
        cellCount = cellCounts(cC);
        mergeOpts.loadFilename  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_aff.mat'];
        mergeOpts.saveFileName  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1'];
        mergeSupervoxels_subgraph(mergeOpts);
        mergeOpts2.loadFilename = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
        mergeOpts2.saveFileName = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2'];
        mergeSupervoxels_subgraph(mergeOpts2);
    end
    
    cellCount = cellCounts(2);
    for bn = 2:length(baseline_noises)
        baseline_noise = baseline_noises(bn);
        mergeOpts.loadFilename  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_aff.mat'];
        mergeOpts.saveFileName  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1'];
        mergeSupervoxels_subgraph(mergeOpts);
        mergeOpts2.loadFilename = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
        mergeOpts2.saveFileName = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2'];
        mergeSupervoxels_subgraph(mergeOpts2);
    end
    
    baseline_noise = baseline_noises(1);
    for pn = 2:length(protein_noises)
        protein_noise = protein_noises(pn);
        mergeOpts.loadFilename  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_aff.mat'];
        mergeOpts.saveFileName  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1'];
        mergeSupervoxels_subgraph(mergeOpts);
        mergeOpts2.loadFilename = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
        mergeOpts2.saveFileName = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2'];
        mergeSupervoxels_subgraph(mergeOpts2);
    end
    
    protein_noise = protein_noises(1);
    for pd = 2:length(protein_densities)
        protein_density = protein_densities(pd);
        mergeOpts.loadFilename  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_aff.mat'];
        mergeOpts.saveFileName  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1'];
        mergeSupervoxels_subgraph(mergeOpts);
        mergeOpts2.loadFilename = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
        mergeOpts2.saveFileName = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2'];
        mergeSupervoxels_subgraph(mergeOpts2);
    end
end

end