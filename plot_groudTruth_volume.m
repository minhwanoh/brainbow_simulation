clear all;
filename = 'sim_9cells_4ch_500bn_0pn_1e-06pd_1.0ef';
load([filename '_denoised.mat'])
load([filename '_6conn_denoised_ws0.01_isplit0.5_20_augmented0.1_merge2_sAff10.mat'])
%load(['simData/sim_' num2str(cellCount) '_ch' num2str(channelCount) '_index.mat'])
bbVol = overallRawVolume;
bbVol(bbVol<0)=0;
for kk = 1:size(bbVol, 4); 
    rawStack = bbVol(:,:,:,kk); 
    rawStack = rawStack - min(rawStack(:)); 
    rawStack = rawStack / max(rawStack(:)); 
    bbVol(:,:,:,kk) = rawStack; 
end; 
clear rawStack;
nonemptyVL = find(~cellfun(@isempty,volumeLabels));
volumeLabels = volumeLabels(nonemptyVL);

seg2=zeros(size(denoised,1),size(denoised,2),size(denoised,3));
for kk=1:numel(volumeLabels)
    seg2(volumeLabels{kk})=kk;
end;
segmentation = seg2;
clusterCount = numel(volumeLabels);
voxelCount                                   = numel(segmentation);

for kk=1:clusterCount
  im1 = zeros(size(segmentation,1), size(segmentation,2), 3);
  tmp2                                       = find(segmentation==kk);
  for mm = 1:3
    tmp1                                     = zeros(size(segmentation));
    tmp1(tmp2)                               = bbVol(tmp2+(mm-1)*voxelCount);
    im_temp = max(tmp1, [], 3);
    im1(:,:,mm) = im_temp;
  end
  %imshow(squeeze(im1),[])
  imwrite(im1, [filename '_neuron' num2str(kk) '.png'])
end

%figure;imshow(squeeze(max(overallRawVolume(:,:,:,1:3),[],3)),[])
vol = squeeze(max(overallRawVolume(:,:,:,1:3),[],3));
imwrite(vol, [filename '_volume.png']);