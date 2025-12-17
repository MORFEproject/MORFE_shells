function PlotBifurcation(Om_HB, Ualpmlitude, StableSystem, Dof, Oms, Ome)
    IterationNum = length(Om_HB);
    hold on
    % bifurcation type
    h_NS  = plot(nan, nan, '^', 'markersize', 10, 'markerfacecolor', [0.4, 0.6, 0.5]); % Neimark-Sacker bifurcation
    h_SN  = plot(nan, nan, 'o', 'markersize', 10, 'markerfacecolor', [0.85, 0.4, 0.2]); % Saddle-node bifurcation
    h_BP  = plot(nan, nan, 's', 'markersize', 10, 'markerfacecolor', [0.6, 0.4, 0.6]); % Branch point bifurcation

    temp = Ualpmlitude(Dof, :); 
    dy_overall = 0.01 * (max(temp(:)) - min(temp(:)));
    dx = 0.01 * (Ome - Oms);
    
    % 
    for i = 1 : IterationNum
        x = Om_HB(i);
        y = Ualpmlitude(Dof, i); 
        markerLabel = ''; 

        switch StableSystem{i}.BifurcationType
            case 'Neimark-Sacker bifurcation'
                plot(x, y, '^', 'markersize', 10, 'markerfacecolor', [0.4, 0.6, 0.5])
                markerLabel = 'NS';
            case 'Saddle-node bifurcation'
                plot(x, y, 'o', 'markersize', 10, 'markerfacecolor', [0.85, 0.4, 0.2])
                markerLabel = 'SN';
            case 'Branch point bifurcation'
                plot(x, y, 's', 'markersize', 10, 'markerfacecolor', [0.6, 0.4, 0.6])
                markerLabel = 'BP';
        end

        if ~isempty(markerLabel)
            for j = 1:length(Dof)
                text(x + dx, y(j) + dy_overall, markerLabel, 'FontSize', 8, 'Color', [0.3, 0.3, 0.3]);
            end
        end
    end

    % 
    lineColor = [0.3, 0.3, 0.7]; 
    h_main_stable = [];
    h_main_unstable = [];
    i = 1;
    while i <= IterationNum
        currentStability = StableSystem{i}.StableIndex; 
        segIdx = i;
        while i+1 <= IterationNum && StableSystem{i+1}.StableIndex == currentStability
            i = i + 1;
            segIdx = [segIdx, i];
        end
        if currentStability == 0
            h = plot(Om_HB(segIdx), Ualpmlitude(Dof, segIdx), '--', 'linewidth', 2, 'color', lineColor);
            if isempty(h_main_unstable)
                h_main_unstable = h;
            end
        else
            h = plot(Om_HB(segIdx), Ualpmlitude(Dof, segIdx), '-', 'linewidth', 2, 'color', lineColor);
            if isempty(h_main_stable)
                h_main_stable = h;
            end
        end
        i = i + 1;
    end

    xlim([Oms Ome])

    h_all = [h_NS(:)', h_SN(:)', h_BP(:)', h_main_unstable(:)', h_main_stable(:)'];

    labels = {'Neimark-Sacker bifurcation', 'Saddle-node bifurcation', ...
              'Branch point bifurcation', 'Unstable system', 'Stable system'};
    legend(h_all, labels, 'Location', 'Best');
end
