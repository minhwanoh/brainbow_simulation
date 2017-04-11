function simulate_par_mergeSupervoxels_subgraph(paramID)

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

mergeOpts.loadFilename  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_aff.mat'];
mergeOpts.saveFileName  = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1'];
mergeSupervoxels_subgraph(mergeOpts);
mergeOpts2.loadFilename = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge1_sAff2.mat'];
mergeOpts2.saveFileName = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2'];
mergeSupervoxels_subgraph(mergeOpts2);

end