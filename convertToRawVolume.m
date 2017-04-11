cellCounts = [5 9 13]; 
for cC = 1:3 
    cellCount = cellCounts(cC);
    for channelCount = 3:6
        clear volumeLabels;
        volfilename012 = ['/Users/min-hwan/Documents/SimExm/output/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_channels_012 .tif'];
        info    = imfinfo(volfilename012);
        thisVol012 = zeros(info(1).Height, info(1).Width, numel(info), 3);
        
        for zz = 1:numel(info);
          temp = imread(volfilename012,'Index',zz);
          thisVol012(:,:,zz,:) = reshape(temp, info(1).Height, info(1).Width, 1, 3);
        end
        overallRawVolume012 = thisVol012/max(thisVol012(:));
        
        if channelCount >= 4
            volfilename345 = ['/Users/min-hwan/Documents/SimExm/output/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_channels_345 .tif'];
            info    = imfinfo(volfilename345);
            thisVol345 = zeros(info(1).Height, info(1).Width, numel(info), 3);

            for zz = 1:numel(info);
              temp = imread(volfilename345,'Index',zz);
              thisVol345(:,:,zz,:) = reshape(temp, info(1).Height, info(1).Width, 1, 3);
            end
            overallRawVolume345 = thisVol345/max(thisVol345(:));
        end
        dim = size(overallRawVolume012);
        overallRawVolume = zeros([dim(1:3), channelCount]);
        overallRawVolume(:,:,:,1:3) = overallRawVolume012;
        if channelCount >= 4
            overallRawVolume(:,:,:,4:channelCount) = overallRawVolume345(:,:,:,1:(channelCount-3));
        end
        
        gtfilename = ['/Users/min-hwan/Documents/SimExm/output/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_gt.tif'];
        info    = imfinfo(gtfilename);
        thisVol = zeros(info(1).Height, info(1).Width, numel(info));
        for zz = 1:numel(info);
          thisVol(:,:,zz) = imread(gtfilename,'Index',zz);
        end
        gt = thisVol;
        k = 1;
        for kk = 1:max(gt(:))
           if sum(gt(:)==kk) > 0
                volumeLabels{k} = (gt==kk);
                k = k + 1;
           end
        end
        colorMatrix = zeros(size(volumeLabels,2), size(overallRawVolume,4));
        save(['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_raw.mat'], 'overallRawVolume', 'volumeLabels','colorMatrix');
    end
end

%imshow(squeeze(max(overallRawVolume(:,:,:,1:3), [], 3)),[]);
