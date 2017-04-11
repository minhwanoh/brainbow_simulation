function convert_raw_volume(cellCount, channelCount, baseline_noise, protein_noise, protein_density, expansion_factor)
    clear volumeLabels;
    volfilename012 = ['/Users/min-hwan/Documents/SimExm/output/' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef/simulation/channels_012.tiff' ];
    info    = imfinfo(volfilename012);
    thisVol012 = zeros(info(1).Height, info(1).Width, numel(info), 3);

    for zz = 1:numel(info)
      temp = imread(volfilename012,'Index',zz);
      thisVol012(:,:,zz,:) = reshape(temp, info(1).Height, info(1).Width, 1, 3);
    end

    for k = 1:3
        tempChannel = thisVol012(:,:,:,k);
        if max(tempChannel(:)) == 0
            display('Empty Channel Detected')
            display([num2str(cellCount) 'cells_' num2str(channelCount) 'channels_' num2str(k) 'ch' ])
        end
    end

    overallRawVolume012 = thisVol012/max(thisVol012(:));

    if channelCount >= 4
        volfilename345 = ['/Users/min-hwan/Documents/SimExm/output/' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef/simulation/channels_345.tiff' ];
        info    = imfinfo(volfilename345);
        thisVol345 = zeros(info(1).Height, info(1).Width, numel(info), 3);

        for zz = 1:numel(info)
          temp = imread(volfilename345,'Index',zz);
          thisVol345(:,:,zz,:) = reshape(temp, info(1).Height, info(1).Width, 1, 3);
        end

        for k = 1:(channelCount-3)
            tempChannel = thisVol345(:,:,:,k);
            if max(tempChannel(:)) == 0
                display('Empty Channel Detected')
                display([num2str(cellCount) 'cells_' num2str(channelCount) 'channels_' num2str(k+3) 'ch' ])
            end
        end

        overallRawVolume345 = thisVol345/max(thisVol345(:));


    end
    dim = size(overallRawVolume012);
    overallRawVolume = zeros([dim(1:3), channelCount]);
    overallRawVolume(:,:,:,1:3) = overallRawVolume012;
    if channelCount >= 4
        overallRawVolume(:,:,:,4:channelCount) = overallRawVolume345(:,:,:,1:(channelCount-3));
    end

    gtfoldername = ['/Users/min-hwan/Documents/SimExm/output/' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef/groundtruth/' ];
    f = dir(gtfoldername);
    for j = 3:length(f)
        ff = dir([f(j).folder '/' f(j).name]);
        if length(ff) == 2
            display('Empty Channel Detected')
            display([num2str(cellCount) 'cells_' num2str(channelCount) 'channels_' f(j).name ])
        end
        for jj = 3:length(ff)
            gtfilename = [ff(jj).folder '/' ff(jj).name];

            info    = imfinfo(gtfilename);
            thisVol = zeros(info(1).Height, info(1).Width, numel(info));
            for zz = 1:numel(info)
              thisVol(:,:,zz) = imread(gtfilename,'Index',zz);
            end
            gt = thisVol;
            k = max(gt(:));
            volumeLabels{k} = (gt==k);
        end
    end

    colorMatrix = zeros(size(volumeLabels,2), size(overallRawVolume,4));
    save(['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sims_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(baseline_noise) 'bn_' num2str(protein_noise) 'pn_' num2str(protein_density) 'pd_' expansion_factor 'ef_raw.mat'], 'overallRawVolume', 'volumeLabels','colorMatrix');
end