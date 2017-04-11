function simulate_par_tissue(paramID)

channelCounts = [3 4 5 6];
cellCounts = [5 9 13];

[kk1 kk2] =ind2sub([numel(channelCounts), numel(cellCounts)], paramID);
paramIndx = [kk1 kk2];
channelCount = channelCounts(paramIndx(1));
cellCount = cellCounts(paramIndx(2));

config = brainbowConfig(channelCount, cellCount);
[overallRawVolume volumeLabels colorMatrix] = brainbowSimulation_3d_raw(config);  
save(['/vega/stats/users/mo2499/bbSimulation/generated_tissues/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' 'generated_tissue.mat'], 'overallRawVolume','volumeLabels','colorMatrix');
dim = size(overallRawVolume);
for i=1:length(volumeLabels)
    tempVol = zeros(dim(1:3));
    tempVol(volumeLabels{i}(:)==1) = i;
    % 32-bit Tiff
    data = uint32(tempVol);
    outputFileName = ['/vega/stats/users/mo2499/bbSimulation/generated_tissues/' num2str(cellCount), 'cells_', num2str(channelCount), 'ch' '/simVol_', num2str(cellCount), 'cells_', num2str(channelCount), 'ch_cellNum'  num2str(i), '.tif'];
    for j=1:dim(3)
        if (j==1)
            t = Tiff(outputFileName,'w');
            tagstruct.ImageLength     = size(data,1);
            tagstruct.ImageWidth      = size(data,2);
            tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample   = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.RowsPerStrip    = 16;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software        = 'MATLAB';
            t.setTag(tagstruct)
            t.write(data(:,:,j));
        else
            t = Tiff(outputFileName,'a');
            tagstruct.ImageLength     = size(data,1);
            tagstruct.ImageWidth      = size(data,2);
            tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample   = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.RowsPerStrip    = 16;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software        = 'MATLAB';
            t.setTag(tagstruct)
            t.write(data(:,:,j))
        end
        t.close();
    end
end

end