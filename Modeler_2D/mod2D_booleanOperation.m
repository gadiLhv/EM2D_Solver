function newPolStruct = mod2D_booleanOperation(pol1,pol2,opString)
% Receives two polygon structures (TBD: Circles, TBD: Curves)
% and performs a boolean operation

polCoors1 = [pol1.x pol1.y];
polCoors2 = [pol2.x pol2.y];

try
[newPol,numberOfLoops] = clipPolygon(polCoors1,...  
                                     polCoors2,...
                                     mod2D_opStrToCode(opString));
catch(lasterror)
  % Catch the latest error message like this:
  fprintf(1,'Error: %s\n',lasterror.message);
end
% TODO: Support circles
% TODO: Suuport closed curves


% Create polygon structure and update number of loops
newPolStruct = mod2D_createPolygonStruct(newPol(1:end-1,1),newPol(1:end-1,2));
newPolStruct.nParts = numberOfLoops;

end