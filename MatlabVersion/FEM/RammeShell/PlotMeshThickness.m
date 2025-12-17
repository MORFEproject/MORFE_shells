function PlotMeshThickness(Mesh, Scale)
    if nargin < 2
        Scale = 1;
    end

    element = Mesh.Element;
    node = Mesh.Node;
    fixnode = Mesh.FixNode;
    forcenode = Mesh.ForceNode;
    freenode = Mesh.FreeNode;
    elementum = Mesh.ElementNum;
    face = [1 2 3 4];
    
    figure;
    hold on;
    for i = 1:elementum
        n1 = element(i,1);
        n2 = element(i,2);
        n3 = element(i,3);
        n4 = element(i,4);
        
        v1 = [node(n1,1), node(n1,2), node(n1,3)];
        v2 = [node(n2,1), node(n2,2), node(n2,3)];
        v3 = [node(n3,1), node(n3,2), node(n3,3)];
        v4 = [node(n4,1), node(n4,2), node(n4,3)];
        
        nVec = cross(v2 - v1, v3 - v1);
        if norm(nVec) == 0
            nVec = [0 0 1]; 
        else
            nVec = nVec / norm(nVec);
        end
        
        t1 = node(n1,4);
        t2 = node(n2,4);
        t3 = node(n3,4);
        t4 = node(n4,4);
        
        v1_top = v1 + 0.5 * Scale * t1 * nVec;
        v1_bot = v1 - 0.5 * Scale * t1 * nVec;
        v2_top = v2 + 0.5 * Scale * t2 * nVec;
        v2_bot = v2 - 0.5 * Scale * t2 * nVec;
        v3_top = v3 + 0.5 * Scale * t3 * nVec;
        v3_bot = v3 - 0.5 * Scale * t3 * nVec;
        v4_top = v4 + 0.5 * Scale * t4 * nVec;
        v4_bot = v4 - 0.5 * Scale * t4 * nVec;
        
        vert_top = [v1_top; v2_top; v3_top; v4_top];
        vert_bot = [v1_bot; v2_bot; v3_bot; v4_bot];
        
        % top-surface
        patch('Faces', face, 'Vertices', vert_top, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 1);
        % bottom-surface
        patch('Faces', face, 'Vertices', vert_bot, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 1);
        
        % 
        vert8 = [v1_top; v2_top; v3_top; v4_top; v1_bot; v2_bot; v3_bot; v4_bot];
        % side surface
        sideFaces = [1 2 6 5;
                     2 3 7 6;
                     3 4 8 7;
                     4 1 5 8];
        patch('Faces', sideFaces, 'Vertices', vert8, 'FaceColor', [0.9290 0.6940 0.1250], 'FaceAlpha', 1);
    end
    
    % 
    for i = 1:length(fixnode)
        idx = fixnode(i);
        x = node(idx,1); y = node(idx,2); z = node(idx,3);
        scatter3(x, y, z, 40, 'red', 'filled');
    end
    for i = 1:length(forcenode)
        idx = forcenode(i);
        x = node(idx,1); y = node(idx,2); z = node(idx,3);
        scatter3(x, y, z, 40, 'yellow', 'filled');
    end
    axis equal;
    axis off;
    view(3);
    hold off;
end
