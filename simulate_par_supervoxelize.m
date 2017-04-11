function simulate_par_supervoxelize(paramID)

superVoxelOpts.brightnessThreshold                                       = 0.1;
superVoxelOpts.spatialDistanceCalculationOpts.upperBound                 = 2;
superVoxelOpts.splitInconsistentSVopts.maxPerimeter                      = 0.5;
superVoxelOpts.splitInconsistentSVopts.connectivity                      = 26;
superVoxelOpts.splitInconsistentSVopts.subdivisionSizeThreshold          = 20;
superVoxelOpts.HMINTH26                                                  = 0.01;

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

superVoxelOpts.dataset = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'];
superVoxelOpts.filePreamble = ['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1'];
cd /vega/stats/users/mo2499/bbSimulation
supervoxelize(superVoxelOpts);

end