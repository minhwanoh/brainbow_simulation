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
        channelCount = channelCounts(ch);
        convert_raw_volume(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)
    end
end

cellCount = cellCounts(2);
for bn = 2:length(baseline_noises)
    baseline_noise = baseline_noises(bn);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch);
        convert_raw_volume(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)
    end
end

baseline_noise = baseline_noises(1);
for pn = 2:length(protein_noises)
    protein_noise = protein_noises(pn);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch);
        convert_raw_volume(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)
    end
end

protein_noise = protein_noises(1);
for pd = 2:length(protein_densities)
    protein_density = protein_densities(pd);
    for ch = 1:length(channelCounts)
        channelCount = channelCounts(ch);
        convert_raw_volume(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)
    end
end