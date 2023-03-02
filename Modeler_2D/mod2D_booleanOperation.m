function newPolStruct = mod2D_booleanOperation(pol1,pol2,opString)
    % Receives two polygon structures (TBD: Circles, TBD: Curves)
    % and performs a boolean operation

    % @item 0: difference @var{inpol} - @var{clippol}
    %
    % @item 1: intersection ("AND") of @var{inpol} and @var{clippol} (= default)
    %
    % @item 2: exclusiveOR ("XOR") of @var{inpol} and @var{clippol}
    %
    % @item 3: union ("OR") of @var{inpol} and @var{clippol}
    % @end itemize
    %

    polCoors1 = [pol1.x pol1.y];
    polCoors2 = [pol2.x pol2.y];

    opCode = mod2D_opStrToCode(opString);
%    try
%    [newPol,numberOfLoops] = clipPolygon(polCoors1,...
%                                         polCoors2,...
%                                         opCode,...
%                                         'clipper');

    [newPol,numberOfLoops] = clipPolygon_clipper(...
                                 polCoors1,...
                                 polCoors2,...
                                 opCode);
    % Check for zero volume polygons

%    catch(lasterror)
%        % Catch the latest error message like this:
%        fprintf(1,'Error: %s\n',lasterror.message);
%    end
    % TODO: Support circles
    % TODO: Suuport closed curves


    % Create polygon structure and update number of loops

    newPolStruct = mod2D_createPolygonStruct(newPol(1:(end-1),1),newPol(1:(end-1),2));
    newPolStruct.nParts = numberOfLoops;

end
