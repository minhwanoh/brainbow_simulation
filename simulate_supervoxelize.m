function simulate_supervoxelize

superVoxelOpts.brightnessThreshold                                       = 0.1;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound                 = 2;
superVoxelOpts.splitInconsistentSVopts.maxPerimeter                      = 0.5;
superVoxelOpts.splitInconsistentSVopts.connectivity                      = 26;
superVoxelOpts.splitInconsistentSVopts.subdivisionSizeThreshold          = 20;
superVoxelOpts.HMINTH26                                                  = 0.01;


baseline_noises = [0, 200, 500];
protein_noises = [0, 0.1, 0.2];
protein_densities = [1e-6, 1e-7, 1e-8];
channelCounts = [3 4 5 6];
cellCounts = [5 9 13];

expansion_factor = '1.0';
baseline_noise = baseline_noises(1);
protein_noise = protein_noises(1);
protein_density = protein_densities(1);

for cC = 1:length(cellCounts)
    cellCount = cellCounts(cC);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        superVoxelOpts.dataset = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'];
        superVoxelOpts.filePreamble = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1'];
        cd /vega/stats/users/mo2499/bbSimulation
        supervoxelize(superVoxelOpts);
    end
end

cellCount = cellCounts(2);
for bn = 2:length(baseline_noises)
    baseline_noise = baseline_noises(bn);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        superVoxelOpts.dataset = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'];
        superVoxelOpts.filePreamble = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1'];
        cd /vega/stats/users/mo2499/bbSimulation
        supervoxelize(superVoxelOpts);
    end
end

baseline_noise = baseline_noises(1);
for pn = 2:length(protein_noises)
    protein_noise = protein_noises(pn);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        superVoxelOpts.dataset = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'];
        superVoxelOpts.filePreamble = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1'];
        cd /vega/stats/users/mo2499/bbSimulation
        supervoxelize(superVoxelOpts);
    end
end

protein_noise = protein_noises(1);
for pd = 2:length(protein_densities)
    protein_density = protein_densities(pd);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        superVoxelOpts.dataset = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'];
        superVoxelOpts.filePreamble = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1'];
        cd /vega/stats/users/mo2499/bbSimulation
        supervoxelize(superVoxelOpts);
    end
end

end

