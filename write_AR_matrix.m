for ii = 1:3
    load(['set_' num2str(ii) '_AR_results.mat']);
    AR{ii} = ARfg_agg;
end

channelCounts = [3 4 5 6];

ar_array = zeros(12*length(AR),3, 4);
for j = 1:4
    for i = 1:length(AR)
        temp = AR{i}(:,:,j);
        for ind = 1:12
            [kk1 kk2] =ind2sub([4, 3], ind);
            ar_array( (ind + (i-1)*12),:,j) = [kk2 temp(kk1,kk2) channelCounts(kk1)];
        end       
    end
end

csvwrite('AR_cellCount.csv', ar_array(:,:,1));
csvwrite('AR_baseline_noise.csv', ar_array(:,:,2))
csvwrite('AR_protein_noise.csv', ar_array(:,:,3))
csvwrite('AR_protein_density.csv', ar_array(:,:,4))
