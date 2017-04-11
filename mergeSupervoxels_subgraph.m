function mergeSupervoxels_subgraph(opts)

thisOpts                                                                 = opts;
load(opts.loadFilename);
opts                                                                     = thisOpts;
if exist('superVoxelCells')
  svCells                                                                = superVoxelCells; clear superVoxelCells;
  svMeans                                                                = superVoxelMeans; clear superVoxelMeans;
  cc                                                                     = numel(svCells);
  [ii_sAff, ~, ss_sAff]                                                  = find(sAff);
  yy                                                                     = ceil( (2*cc-1 - sqrt((2*cc-1)^2-8*ii_sAff))/2 );
  xx                                                                     = ii_sAff - cc*yy + cc + yy.*(yy+1)/2;
  square_sAff                                                            = sparse(xx, yy, ss_sAff, cc, cc);
  square_sAff                                                            = square_sAff + transpose(square_sAff);
  clear ii_sAff; clear ss_sAff; clear xx; clear yy; clear sAff;
end
voxelCounts                                                              = cellfun(@numel, svCells);
voxelCount                                                               = prod(stackSize);
allTriplets                                                              = nchoosek(1:size(svMeans, 2), 3);
allColors                                                                = zeros(size(svMeans, 1), 3*size(allTriplets, 1));
for kk = 1:size(allTriplets, 1)
  allColors(:, 3*kk-2:3*kk)                                              = rgb2luv(svMeans(:, allTriplets(kk, :))')';
end
[coeff,score,latent]                                                     = pca(allColors);
svMeansLUV                                                               = score(:, 1:size(svMeans, 2));
tic; [vdpgmAssignments clusterCount]                                     = vdpgmEstimate(svMeansLUV); disp(clusterCount); toc;
binsaff                                                                  = (square_sAff > 1/(sqrt(3)+1e-4));
[rows, cols]                                                             = find(binsaff);
differentColors                                                          = find(sum((svMeansLUV(rows, :) - svMeansLUV(cols, :)).^2, 2) > opts.kmeansMerging.cutSubgraphsMaxLuvColorDistance^2);
rows(differentColors)                                                    = [];
cols(differentColors)                                                    = [];
binsaff                                                                  = sparse(rows, cols, 1, size(binsaff, 1), size(binsaff, 2));
subgraphNodes                                                            = randomlyPartitionBinaryGraph(binsaff, opts.kmeansMerging.maxSubgraphSize);
svCount                                                                  = cell(1, numel(subgraphNodes));
stackSize                                                                = stackSize; % for parfor
parfor kk = 1:numel(subgraphNodes)
  optskk = opts; oldCountkk = numel(subgraphNodes{kk}); svCount{kk} = [];
  nkk = subgraphNodes{kk}; sAffkk = square_sAff(nkk, nkk); svMeanskk = svMeans(nkk, :); svCellskk = svCells(nkk); vCkk = voxelCounts(nkk); svCMinskk = svColorMins(nkk, :); svCMaxskk = svColorMaxs(nkk, :); vdpgmkk = vdpgmAssignments(nkk);
  optskk.kmeansMerging.clusterCount                                      = optskk.kmeansMerging.overSamplingFactor * numel(unique(vdpgmAssignments(nkk)));
  oI                                                                     = 1:numel(nkk);
  while true
    cd /vega/stats/users/us2157/bb/bbSimulation;
    optskk.mergeSingleNeighborSuperVoxels.maxSizeForSingleNeighborSVs    = quantile(vCkk, 0.5);
    % HEURISTICS TO MERGE SUPERVOXELS
    [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, ~]    = demixSupervoxels(                sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, opts.demix, oI);
    [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI]       = mergeWRTneighborsAndOrientations(sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, opts.mergeWRTnAo.normFlag, stackSize, opts.mergeWRTnAo);
    [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI]       = mergeSingleNeighborSuperVoxels(  sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, optskk.mergeSingleNeighborSuperVoxels);
    [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI]       = mergeCloseNeighborhoods(         sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, opts.mergeCloseNeighborhoods);
    [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI]       = mergeSmallSuperVoxels(           sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, opts.mergeSmallSuperVoxels);
    if optskk.kmeansMerging.clusterCount<numel(svCellskk)
      [sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI]     = mergeWithKmeans(                 sAffkk, svMeanskk, svCellskk, vCkk, svCMinskk, svCMaxskk, oI, optskk.kmeansMerging);
    end
    svCount{kk}(end+1) = numel(svCellskk);
    if svCount{kk}(end)>oldCountkk-5; subgraph_svColorMins{kk} = svCMinskk; subgraph_svColorMaxs{kk} = svCMaxskk; subgraph_svCells{kk} = svCellskk; subgraph_indices{kk} = oI; break; end;
    oldCountkk = svCount{kk}(end);
  end
end
[svMeans, svCells, voxelCounts,  svColorMins, svColorMaxs, idxTransform] = mergeSubgraphs(svMeans, voxelCounts, subgraph_svColorMins, subgraph_svColorMaxs, subgraph_svCells, subgraph_indices, subgraphNodes);
disp(numel(svCells))
% SPATIAL AFFINITY CALCULATION
square_sAff                                                              = calculate_square_sAff(svCells, stackSize, opts.spatialDistanceCalculation, opts.zAnisotropy);
thisFN                                                                   = [opts.saveFileName '_sAff' num2str(opts.spatialDistanceCalculation.upperBound) '.mat'];
save(thisFN, 'superVoxelOpts', 'opts', 'svCells', 'svMeans', 'voxelCounts', 'stackSize', 'square_sAff', 'boundaryVoxels', 'svCount', 'svColorMins', 'svColorMaxs', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, C] = mergeSubgraphs(svMeans, voxelCounts, subgraph_svColorMins, subgraph_svColorMaxs, subgraph_svCells, subgraph_indices, nodePartition)
% nodePartition: cell variable holding the supervoxels given to individual subsets of the node partition
% subgraph_indices: cell variable holding the index variables for individual subsets of the node partition (returned by mergeConnectedComponentsOfSubgraph)
% subgraph_svCells: cell variable holding the supervoxel cells for individual subsets of the node partition (returned by mergeConnectedComponentsOfSubgraph)
oldsvCount                                                       = size(svMeans, 1);
newsvCount                                                       = sum(cellfun(@numel, subgraph_svCells));
svCells                                                          = cell(1, newsvCount);
svColorMins                                                      = zeros(newsvCount, size(svMeans, 2));
svColorMaxs                                                      = zeros(newsvCount, size(svMeans, 2));
C                                                                = zeros(1, oldsvCount);
offset                                                           = 0;
voxelCounts                                                      = voxelCounts(:);
for kk = 1:numel(subgraph_indices)
  C(nodePartition{kk})                                           = offset + subgraph_indices{kk};
  svCells(offset+1:offset+numel(subgraph_svCells{kk}))           = subgraph_svCells{kk};
  svColorMins(offset+1:offset+numel(subgraph_svCells{kk}), :)    = subgraph_svColorMins{kk};
  svColorMaxs(offset+1:offset+numel(subgraph_svCells{kk}), :)    = subgraph_svColorMaxs{kk};
  offset                                                         = offset + numel(subgraph_svCells{kk});
end
newsvMeans                                                       = zeros(newsvCount, size(svMeans, 2));
for kk = 1:newsvCount
  thisGroup                                                      = find(C==kk);
  newsvMeans(kk, :)                                              = voxelCounts(thisGroup)' * svMeans(thisGroup, :) / sum(voxelCounts(thisGroup));
end
svMeans                                                          = newsvMeans;
voxelCounts                                                      = cellfun(@numel, svCells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function subgraphNodes = randomlyPartitionBinaryGraph(binsaff, maxSubgraphSize)
[S, C]                   = graphconncomp(binsaff);
subgraphNodes            = cell(0);
for kk = 1:S
  theseNodes             = find(C==kk);
  if numel(theseNodes) <= maxSubgraphSize
    subgraphNodes{end+1} = theseNodes;
  end
end
availableNodes           = setdiff(1:size(binsaff, 1), cell2mat(subgraphNodes));
while ~isempty(availableNodes)%true
  thisSeed               = randi(numel(availableNodes));
  tmp                    = growSubgraphOnBinaryGraph(binsaff(availableNodes, availableNodes), thisSeed, maxSubgraphSize);
  subgraphNodes{end+1}   = availableNodes(tmp);
  availableNodes         = availableNodes(setdiff(1:numel(availableNodes), tmp, 'stable'));
  if numel(availableNodes)<maxSubgraphSize
    subgraphNodes{end+1} = availableNodes;
    break;
  end
end

while true
  allCounts                 = cellfun(@numel, subgraphNodes);
  [allCounts, idx]          = sort(allCounts);
  last                      = find(cumsum(allCounts)>maxSubgraphSize, 1) - 1;
  if isempty(last)
    newSubgraphNodes{1}     = cell2mat(subgraphNodes);
    break;
  end
  if last > 1
    newSubgraphNodes        = cell(1, numel(subgraphNodes)-last+1);
    newSubgraphNodes{1}     = cell2mat(subgraphNodes(idx(1:last)));
    newSubgraphNodes(2:end) = subgraphNodes(idx(last+1:end));
    subgraphNodes           = newSubgraphNodes;
  else
    break;
  end
end

disp(['# subgraphs: ', num2str(numel(subgraphNodes))]);
%disp(cellfun(@numel, subgraphNodes));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeSubset = growSubgraphOnBinaryGraph(binsaff, node, maxSubgraphSize)
nodeSubset      = node;
newlyAdded      = node;
oldCount        = 1;
while true
  newlyAdded    = setdiff(find(any(binsaff(newlyAdded, :), 1)), nodeSubset, 'stable');
  nodeSubset    = [nodeSubset newlyAdded];
  count         = numel(nodeSubset);
  if count==oldCount | count>=maxSubgraphSize
    break;
  else
    oldCount    = count;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function square_sAff = calculate_square_sAff(svCells, stackSize, spatialDistanceCalculation, zAnisotropy)
cc                                                       = numel(svCells);
boundaryVoxelsSub                                        = cell(1, cc);
parfor kk = 1:cc
  [xx,yy,zz]                                             = ind2sub(stackSize, svCells{kk});
  xSub                                                   = min(xx)-2;
  ySub                                                   = min(yy)-2;
  zSub                                                   = min(zz)-2;
  xx                                                     = xx-xSub;
  yy                                                     = yy-ySub;
  zz                                                     = zz-zSub;
  maxxx                                                  = max(xx);
  maxyy                                                  = max(yy);
  maxzz                                                  = max(zz);
  tmp                                                    = false(maxxx+1, maxyy+1, maxzz+1);
  reducedIndices                                         = sub2ind([maxxx+1, maxyy+1, maxzz+1], xx, yy, zz);
  tmp(reducedIndices)                                    = true;
  localBoundaryVoxels                                    = find(tmp & ~imerode(tmp, ones(3,3,3)));
  if ~isempty(localBoundaryVoxels)
    [xx,yy,zz]                                           = ind2sub(size(tmp), localBoundaryVoxels);
  end
  xx                                                     = xx+xSub;
  yy                                                     = yy+ySub;
  zz                                                     = zz+zSub;
  boundaryVoxelsSub{kk}                                  = [xx,yy,zz*zAnisotropy];
end
sAff                                                     = calculate_sAff(cc, boundaryVoxelsSub, spatialDistanceCalculation);
[ii_sAff, ~, ss_sAff]                                    = find(sAff);
yy                                                       = ceil( (2*cc-1 - sqrt((2*cc-1)^2-8*ii_sAff))/2 );
xx                                                       = ii_sAff - cc*yy + cc + yy.*(yy+1)/2;
square_sAff                                              = sparse(xx, yy, ss_sAff, cc, cc);
square_sAff                                              = square_sAff + transpose(square_sAff);
