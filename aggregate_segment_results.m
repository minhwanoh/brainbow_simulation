function aggregate_segment_results

baseline_noises = [0, 200, 500];
protein_noises = [0, 0.1, 0.2];
protein_densities = [1e-6, 1e-7, 1e-8];
channelCounts = [3 4 5 6];
cellCounts = [5 9 13];

expansion_factor = '1.0';

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
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat']);
        new_row = [expressionDensity, ARfg, channelCount]
        AR_mat_cellCount = [AR_mat_cellCount; new_row]; 
    end

    % baseline_noise experiment
    cellCount = cellCounts(2);
    for bn = 1:length(baseline_noises)
        baseline_noise = baseline_noises(bn);
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat']);
        new_row = [baseline_noise, ARfg, channelCount]
        AR_mat_baseline_noise = [AR_mat_baseline_noise; new_row];
    end
    
    % protein_noise experiment
    baseline_noise = baseline_noises(1);
    for pn = 1:length(protein_noises)
        protein_noise = protein_noises(pn);
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat']);
        new_row = [protein_noise, ARfg, channelCount]
        AR_mat_protein_noise = [AR_mat_protein_noise; new_row];
    end
    
    % protein density experiment
    protein_noise = protein_noises(1);
    for pd = 1:length(protein_densities)
        protein_density = protein_densities(pd);
        load(['/vega/stats/users/mo2499/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_index.mat']);
        new_row = [protein_density, ARfg, channelCount]
        AR_mat_protein_density = [AR_mat_protein_density; new_row];
    end
        
end

t = datetime('now');
DateString = datestr(t);
save(['/vega/stats/users/mo2499/bbSimulation/simData/AR_results_' DateString '.mat'], 'AR_mat_cellCount', 'AR_mat_baseline_noise', 'AR_mat_protein_noise', 'AR_mat_protein_density');
end