function [square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex] = mergeWithKmeans(square_sAff, svMeans, svCells, voxelCounts, svColorMins, svColorMaxs, origIndex, opts)

colorTriplets                                                   = nchoosek(1:size(svMeans, 2), 3);
allLUV                                                          = zeros(size(svMeans, 1), size(colorTriplets, 1)*3);
for kk = 1:size(colorTriplets, 1)
  allLUV(:, (kk-1)*3+1:kk*3)                                    = rgb2luv(svMeans(:, colorTriplets(kk, :))')';
end
[coeff,score,latent]                                            = pca(allLUV);
colorsForClustering                                             = score(:, 1:min(size(svMeans, 2), size(score, 2)));
index                                                           = kmeans(colorsForClustering, opts.clusterCount, 'MaxIter', 1000, 'Replicates', 100);
binsaff_6n                                                      = (square_sAff>1-eps);
C                                                               = zeros(size(index));
for kk = 1:max(index)
  theseSVs                                                      = find(index==kk);
  thisLUV                                                       = colorsForClustering(theseSVs, :);
  bincaff                                                       = (squareform(pdist(thisLUV))<opts.maxColorDistance);
  if isempty(bincaff)
    bincaff                                                     = true;
  end
  [thisS, thisC]                                                = graphconncomp(binsaff_6n(theseSVs, theseSVs) & bincaff, 'Weak', true);
  C(theseSVs)                                                   = max(C) + thisC;
end
C(index==0)                                                     = max(C)+1:max(C)+nnz(index==0);
S                                                               = max(C);

newsvCells              = cell(1, S);
newsvMeans              = zeros(S, size(svMeans, 2));
newsvColorMins          = zeros(S, size(svMeans, 2));
newsvColorMaxs          = zeros(S, size(svMeans, 2));

voxelCounts             = voxelCounts(:);

for kk = 1:S
  thisConnComp          = find(C==kk);
  newsvCells{kk}        = cell2mat(svCells(thisConnComp)');
  newsvMeans(kk, :)     = voxelCounts(thisConnComp)' * svMeans(thisConnComp, :) / sum(voxelCounts(thisConnComp));
  newsvColorMins(kk, :) = min(svColorMins(thisConnComp, :), [], 1);
  newsvColorMaxs(kk, :) = max(svColorMaxs(thisConnComp, :), [], 1);
end
svCells                 = newsvCells;
svMeans                 = newsvMeans;
svColorMins             = newsvColorMins;
svColorMaxs             = newsvColorMaxs;
voxelCounts             = cellfun(@numel, newsvCells);

[row, col, val]         = find(square_sAff);
row                     = C(row);
col                     = C(col);
upper                   = find(row<=col);
row(upper)              = [];
col(upper)              = [];
val(upper)              = [];

ff                      = max(row)+1;
id                      = row + col*ff';
[uniqueid, ia, ~]       = unique(id);
newRows                 = row(ia);
newCols                 = col(ia);
newVals                 = zeros(size(newRows));
parfor kk = 1:numel(uniqueid)
  newVals(kk)           = max(val(id==uniqueid(kk)));
end
square_sAff             = sparse(newRows, newCols, newVals, S, S);
square_sAff             = square_sAff + transpose(square_sAff);

%origIndex               = C(origIndex);
origIndex(origIndex~=0) = C(origIndex(origIndex~=0));

