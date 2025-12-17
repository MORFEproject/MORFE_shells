function PlotModalShape(Mesh, Phi, ComputeMode, Start, End)
    % Calculate the number of modes to plot
    numModes = End - Start + 1;
    
    % Calculate the optimal grid layout (as square as possible)
    nRows = floor(sqrt(numModes));
    nCols = ceil(numModes / nRows);
    
    % If the grid still isn't square enough, try to make it more square
    if nCols > nRows + 1
        nRows = nRows + 1;
        nCols = ceil(numModes / nRows);
    end
    
    % Create a single figure
    figure('Name', 'Modal Shapes', 'Position', [100, 100, 1000, 800]);
    
    for i = Start:End
        % Create subplot for this mode
        subplot(nRows, nCols, i - Start + 1);
        
        % Process mode shape
        PhiFull = zeros(Mesh.Sdof, ComputeMode);
        PhiFull(Mesh.Freedof,:) = Phi;
        PhiSelectFull = PhiFull(:, i);
        PhiSelectFull = reshape(PhiSelectFull, 6, []);
        
        % Plot the deformation for this mode
        PlotDeformation(Mesh, PhiSelectFull', 1);
        
        % Add a title to identify the mode number
        title(['Mode ', num2str(i)]);
    end
    
    % Adjust subplot spacing
    sgtitle('Modal Shapes');
    set(gcf, 'Color', 'w');
    tight_subplot = false;
    if exist('tight_subplot', 'file') == 2 || tight_subplot
        % If tight_subplot function is available, use it for better spacing
        % This is optional and depends on having the tight_subplot function
        h = tight_subplot(nRows, nCols, [.03 .03], [.1 .05], [.05 .05]);
        for i = 1:numModes
            set(gcf, 'CurrentAxes', h(i));
            PhiFull = zeros(Mesh.Sdof, ComputeMode);
            PhiFull(Mesh.Freedof,:) = Phi;
            PhiSelectFull = PhiFull(:, Start+i-1);
            PhiSelectFull = reshape(PhiSelectFull, 5, []);
            PlotDeformation(Mesh, PhiSelectFull', 1);
            title(['Mode ', num2str(Start+i-1)]);
        end
    end
end