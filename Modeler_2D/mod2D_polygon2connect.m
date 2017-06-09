function [x,y,conn] = mod2D_polygon2connect(pol)
% [x,y,conn] = mod2D_polygon2connect(pol)
%
% Converts between a polygon (multiple loops supported) and a connectivity list
% Inputs:
% pol - Polygon structure.
% 
% Outputs:
% x,y - [N,1] coordinates, where N is the number of nodes of all the loops
% conn - [Nc,2] connectivity between nodes.



x = [];
y = [];

% Add NaNs so it will be easier to count the number of loops
pol.x = [pol.x ; NaN];
pol.y = [pol.y ; NaN];

% Initialize connectivity list
conn = [];

% Node counter to know when to begin the connectivity list
cNodes = 0;

for loopIdx = 1:pol.nParts
  nanIdxs = find(isnan(pol.x));
  
  % Snip out current loop
  cx = pol.x(1:nanIdxs(1)-1);
  cy = pol.y(1:nanIdxs(1)-1);
  
  % Check if first and last line are the same
  if (cx(1) == cx(end)) && (cy(1) == cy(end))
    cx(end) = [];
    cy(end) = [];
  end
  
  % Create temporary connectivity list
  cConn = [(1:numel(cx)).' [(2:numel(cx)).' ; 1]];
  
  % Advance all nodes in this loop's connectivity list to correct position
  conn = [conn ; cConn+cNodes];
  
  % Advance node counter
  cNodes = cNodes + numel(cx);
  
  x = [x ; cx];
  y = [y ; cy];
  pol.x(1:nanIdxs) = [];
  pol.y(1:nanIdxs) = [];
end
