function addThicknessColumn(inputFile, outputFile, tmin, tmax)
    data = textread(inputFile);  % data: [nodeID, x, y, z]

    % -------------------------------------------------------------------------
    nodeID = data(:,1);
    xCoord = data(:,2);
    yCoord = data(:,3);
    zCoord = data(:,4);

    % -------------------------------------------------------------------------
    xMin = min(xCoord);
    xMax = max(xCoord);

    % -------------------------------------------------------------------------
    %   t(x) = tmin + (tmax - tmin)*( (x - xMin)/(xMax - xMin) )
    thickness = tmin + (tmax - tmin) .* ( (xCoord - xMin) ./ (xMax - xMin) );

    % [nodeID, x, y, z, thickness]
    % -------------------------------------------------------------------------
    outData = [nodeID, xCoord, yCoord, zCoord, thickness];

    dlmwrite(outputFile, outData, 'delimiter', ' ', 'precision', 15);
    fprintf('New node document with thickness: %s\n', outputFile);
end
