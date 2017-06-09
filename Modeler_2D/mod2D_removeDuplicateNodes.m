function [node,connectivity,faceAssignment] = mod2D_removeDuplicateNodes(node,connectivity,faceAssignment,uniqueDistTH)
  
% Get unique node list
%[nodeUnique,~,unique2origIdx] = unique(node,'rows');
[nodeUnique,~,unique2origIdx] = mod2D_uniqueNodeByDist(node,uniqueDistTH);
% nodeUnique - Unique nodes, explicit
% uniqueIdxs - Extracts the unique rows out of 'node'
%  => nodeUnique = node(uniqueIdx,:)
% unique2origIdxs - Rebuilds the 'node' list out of the unique node list 
%  => node = nodeUnique(unique2origIdx,:)

% If necessary, remove duplicate nodes
if(numel(nodeUnique) ~= numel(node))
  % Keep list of original node indexes
  origIdx = 1:size(node,1);
  
  % Build legend between unique indices and original list
  uniqueIdx = 1:size(nodeUnique,1);
  uniqueIdxLegend = uniqueIdx(unique2origIdx).';
  
  % Each index in the CONNECTIVITY list is compared to original index list
  origDims = size(connectivity);
  connectivity = connectivity(:);
  binMapping = bsxfun(@eq,connectivity,origIdx(:).');
  
  % Each column is multiplied by it's unique counterpart
  replacementIdx = bsxfun(@times,binMapping,uniqueIdxLegend(:).');
  [whichConn,whichLegend] = find(binMapping);
  
  % And the replacement, unique, index, is placed in it's place.
  connectivity(whichConn) = replacementIdx(whichConn + (whichLegend-1)*size(binMapping,1));
  
  % Back to the original dimensions
  connectivity = reshape(connectivity,origDims);
  
  % Store to node list
  node = nodeUnique;
end


% Due to distance threshold, it is possible that there are now self-pointing 
% edges. First, they need to be removed from the edge list. Then, from the 
% respective faces.
zeroEdgeIdx = find(connectivity(:,1) == connectivity(:,2));
if ~isempty(zeroEdgeIdx)
  % Remove zero length edges
  connectivity(zeroEdgeIdx,:) = [];
  faceAssignment(zeroEdgeIdx,:) = [];
end


end