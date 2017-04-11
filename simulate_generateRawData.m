function simulate_generateRawData

cellCounts = [2 5 9 13];
for cC = 1:4
  cellCount = cellCounts(cC);
    for channelCount = 3:5
			 randomWalkSD = 0.1;
%      for randomWalkSD = 0.01:0.01:0.04
        for normalSD = 0.05:0.05:0.1
	  config=brainbowConfig; config.cellsToUse = 1:cellCount; config.channelCount = channelCount; config.colors.randomWalkSD = randomWalkSD; config.normalSD = normalSD;
          tic; [overallRawVolume volumeLabels colorMatrix] = brainbowSimulation_3d_raw(config); toc;
%        origZproj = squeeze(max(overallRawVolume,[],3)); origYproj = squeeze(max(overallRawVolume,[],2)); origXproj = squeeze(max(overallRawVolume,[],1));
%        overallRawVolume = overallRawVolume(301:500,301:500,26:125, :); for kk = 1:numel(volumeLabels); volumeLabels{kk} = volumeLabels{kk}(301:500,301:500,26:125); end;
          save(['/vega/stats/users/us2157/bb/bbSimulation/simData/sim_' num2str(cellCount) 'cells_' num2str(channelCount) 'ch_' num2str(randomWalkSD) 'brownian_' num2str(normalSD) 'normalSD_6conn_raw.mat'], 'overallRawVolume','volumeLabels','colorMatrix');
      end
%    end
  end
end
