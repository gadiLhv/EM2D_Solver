function binZeroVolume = isZeroVolume(X,Y);
    % First, check that the first and last points are not the same
    if (X(1) == X(end)) && (Y(1) == Y(end))
        X(end) = [];
        Y(end) = [];
    end
    
    shiftOne = @(v) [v(end) ; v(1:end-1)];
    
    X1 = X(:);
    Y1 = Y(:);
    X2 = shiftOne(X1);
    Y2 = shiftOne(Y1);
    X3 = shiftOne(X2);
    Y3 = shiftOne(Y2);
    
    vol = 0.5*(...
            X1.*(Y2 - Y3) + ...
            X2.*(Y3 - Y1) + ...
            X3.*(Y1 - Y2));
            
    % If any of the triangles is zero volume, then bye bye
    binZeroVolume = sum(vol < 1e-15) ~= 0;
end
