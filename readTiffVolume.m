cellCounts = [5 9 13]; % [2 5 9 13];
for cC = 1:3 % 4
    clear volumeLabels;
    cellCount = cellCounts(cC);
    volfilename = ['/simulatedTiffs/Composite_' num2str(cellCount) 'cells.tif'];
    info    = imfinfo(volfilename);
    thisVol = zeros(info(1).Height, info(1).Width, numel(info));
    for zz = 1:numel(info);
      thisVol(:,:,zz) = imread(volfilename,'Index',zz);
    end
    bbVol = zeros(size(thisVol,1),size(thisVol,2),size(thisVol,3)/4, 4); %%%% 4 CHANNELS
    for kk=1:size(bbVol,3); 
        bbVol(:,:,kk,1)=thisVol(:,:,4*(kk-1)+1); 
        bbVol(:,:,kk,2)=thisVol(:,:,4*(kk-1)+2); 
        bbVol(:,:,kk,3)=thisVol(:,:,4*(kk-1)+3); 
        bbVol(:,:,kk,4)=thisVol(:,:,4*kk); 
    end;
    overallRawVolume = bbVol/max(bbVol(:));

    gtfilename = ['/simulatedTiffs/gt_' num2str(cellCount) 'cells.tif'];
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
    save(['/Users/min-hwan/Documents/Brainbow/bbSimulation/simData/sim_' num2str(cellCount) 'cells_ch4_raw.mat'], 'overallRawVolume', 'volumeLabels','colorMatrix');

end