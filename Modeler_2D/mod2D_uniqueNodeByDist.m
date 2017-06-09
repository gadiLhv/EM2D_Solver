function [uNode,uIdx,u2oIdx] = mod2D_uniqueNodeByDist(node,distTH)


% Create a distance map between all nodes
distMap = permute(node,[1 3 2]);
distMap = bsxfun(@minus,distMap,permute(distMap,[2 1 3]));
distMap = sqrt(sum(distMap.^2,3));


% Set all upper matrix to infinite distances, to reduce the need to:
% a. Remove the node itself (it is definitely always the closest).
% b. Remove the duplicates.
distMap(triu(logical(ones(size(node,1))))) = inf;

% Find within the possibilities all the nodes that are close enough (or
% coinsiding) with the original nodes.
binSameNode = distMap < distTH;

% Find all duplicates.
[duplicateIdx,origIdx] = find(binSameNode);

% Repeat this until no non-unique indices exist
while(1)
  % Check for matches between the duplicates and the originals (possible in case
  % node is duplicated more than once).
  [whichDup,whichOrig] = find(bsxfun(@eq,duplicateIdx,origIdx.'));

  if isempty(whichDup)
    break;
  end
  % Replace the non-unique "original" indices
  origIdx(whichOrig) = origIdx(whichDup);
  
end

% Remove duplicate nodes
uNode = node;
uNode(duplicateIdx,:) = [];

% Generate the index list that creates uNode = node(uIdx)
uIdx = (1:size(node,1)).';
uIdx(duplicateIdx) = [];

% Generate the legend that allows node = uNnode(u2oIdx):
% These two lists are the same length. One denotes the 
% positions of the unique indices in the original node list.
% The one below represents the indices in the unique list.
uRefIdx = (1:numel(uIdx)).';

u2oIdx = zeros([size(node,1) 1]);
% In order rebuild the original node list, the non-unique
% (duplicate) places need to be filled with the duplicate 
% values.

% 1. Fill with known, unique, values
u2oIdx(uIdx) = uRefIdx;

% 2. Create the full unique list. Here, in the non-unique 
% places, the indexes of the unique indices are placed.
u2oIdx(duplicateIdx) = u2oIdx(origIdx);

%%% Verify results (debug)
%%debugDist = sqrt(sum((uNode(u2oIdx,:) - node).^2,2));
%%plot(debugDist);

end
