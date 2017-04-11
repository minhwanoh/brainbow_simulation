function simulate_denoise
parpool('local',16) 

cd /vega/stats/users/us2157/bb/BM4D_v3p2/

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
    parfor ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_raw.mat']);
        denoised = overallRawVolume;
        sigma = 1/8;
        for kk=1:channelCount; 
            [tmp, ~] = bm4d(squeeze(denoised(:,:,:,kk)), 'Gauss', sigma); 
            denoised(:,:,:,kk) = tmp; 
        end;
        save(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'overallRawVolume','denoised','volumeLabels','colorMatrix');
    end
end

cellCount = cellCounts(2);
for bn = 2:length(baseline_noises)
    baseline_noise = baseline_noises(bn);
    parfor ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_raw.mat']);
        denoised = overallRawVolume;
        sigma = 1/8;
        for kk=1:channelCount; 
            [tmp, ~] = bm4d(squeeze(denoised(:,:,:,kk)), 'Gauss', sigma); 
            denoised(:,:,:,kk) = tmp; 
        end;
        save(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'overallRawVolume','denoised','volumeLabels','colorMatrix');
    end
end

baseline_noise = baseline_noises(1);
for pn = 2:length(protein_noises)
    protein_noise = protein_noises(pn);
    parfor ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_raw.mat']);
        denoised = overallRawVolume;
        sigma = 1/8;
        for kk=1:channelCount; 
            [tmp, ~] = bm4d(squeeze(denoised(:,:,:,kk)), 'Gauss', sigma); 
            denoised(:,:,:,kk) = tmp; 
        end;
        save(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'overallRawVolume','denoised','volumeLabels','colorMatrix');
    end
end

protein_noise = protein_noises(1);
for pd = 2:length(protein_densities)
    protein_density = protein_densities(pd);
    parfor ch = 1:length(channelCounts)
        channelCount = channelCounts(ch)
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_raw.mat']);
        denoised = overallRawVolume;
        sigma = 1/8;
        for kk=1:channelCount; 
            [tmp, ~] = bm4d(squeeze(denoised(:,:,:,kk)), 'Gauss', sigma); 
            denoised(:,:,:,kk) = tmp; 
        end;
        save(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_denoised.mat'], 'overallRawVolume','denoised','volumeLabels','colorMatrix');
    end
end

end