function arProcess

unprocessed = true(1, numel(voxels));
values      = zeros(channelCount, numel(voxels));

while any(unprocessed)
  thisPos   = find(unprocessed, 1);
  thisVoxel = voxels(thisPos);
  [xx yy zz] = ind2sub(stackSize, thisVoxel);
  neighbors = repmat([xx yy zz], 27, 1) + meshgrid([-1 0 1], [-1 0 1], [-1 0 1]);
  neighbors(14, :) = [];
  invalid   = find(sum(neighbors<1,2) | neighbors(:,1)>stackSize(1) | neighbors(:,2)>stackSize(2) | neighbors(:,3)>stackSize(3));
  neighbors(invalid, :) = [];
  neighbors = sub2ind(stackSize, neighbors);
  processedNeighbors = neighbors(~unprocessed(neighbors));
  values(:, thisVoxel) = mean(values(:, processedNeighbors)) + options.sd * randn(channelCount, 1);
  unprocessed(thisPos) = false;
end
