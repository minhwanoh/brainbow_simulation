function config = brainbowConfig(channelCount, cellCount)

config.cellsToUse              = randsample(20,cellCount)'; % 1:9; % 1:20; %[4 8 10 12 13 16 17 18 19];
config.ANGLESTEP               = 90;
config.SHIFTSTEP               = 20;
config.ZSHIFTSTEP              = 5;
config.xSize                   = 400; % 800;
config.ySize                   = 400; % 800;
config.zSize                   = 150; % 150;
config.overflow                = 20; % 5;
config.channelCount            = channelCount;
config.colors.preAssignedRatio = 0.1;
config.colors.randomWalkSD     = 0.03;
config.normalSD                = 0.05;
config.maxConnCompPerCell      = 3;
