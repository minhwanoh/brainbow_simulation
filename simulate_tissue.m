cellCounts = [5 9 13];
for cC = 1:3
    cellCount = cellCounts(cC);
    for channelCount = 3:6
        config = brainbowConfig(channelCount, cellCount);
        [overallRawVolume volumeLabels colorMatrix] = brainbowSimulation_3d_raw(config);  
        save(['/vega/stats/users/mo2499/bbSimulation/generated_tissues/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' 'generated_tissue.mat'], 'overallRawVolume','volumeLabels','colorMatrix');
        dim = size(overallRawVolume);
        tempVol = zeros(dim(1:3));
        for i=1:length(volumeLabels);
            tempVol(volumeLabels{i}(:)==1) = i;
        end
        % 32-bit Tiff
        data = uint32(tempVol);
        for i=1:dim(3);
            outputFileName = ['/vega/stats/users/mo2499/bbSimulation/generated_tissues/' num2str(cellCount), 'cells_', num2str(channelCount), 'ch' '/simVol_', num2str(cellCount), 'cells_', num2str(channelCount), 'ch_'  num2str(100+i), '.tif'];
            % This is a direct interface to libtiff
            t = Tiff(outputFileName,'w');
            % Setup tags
            % Lots of info here:
            % http://www.mathworks.com/help/matlab/ref/tiffclass.html
            tagstruct.ImageLength     = size(data,1);
            tagstruct.ImageWidth      = size(data,2);
            tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample   = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.RowsPerStrip    = 16;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software        = 'MATLAB';
            t.setTag(tagstruct)
            t.write(data(:,:,i));
            t.close();
        end
    end
end