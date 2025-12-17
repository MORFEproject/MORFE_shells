function PlotDeformation(Mesh, Up, Scale)
    element = Mesh.Element;
    node = Mesh.Node;
    elementnum = Mesh.ElementNum;
    face = [1 2 3 4];
    % figure
    % title('Deformation')
    hold on 
    for i = 1 : elementnum
        Ue = Up(element(i,:),:);
        Ue_plane = Ue(:, [1:3]);
        we = Ue(:,3);
        wreduce = we(face);
        NodeLocal = node(element(i,:),1:3) + Scale * Ue_plane;
        NodeLocal2 = node(element(i,:),1:3);
        vert = NodeLocal(face,:);
        vert2 = NodeLocal2(face,:);
        patch(vert2(:,1),vert2(:,2),vert2(:,3),[0 0.4470 0.7410],'FaceAlpha',0.2);
        patch(vert(:,1),vert(:,2),vert(:,3),wreduce,'FaceColor','interp');
    end
    axis equal; axis off; box off;
    colormap;
    
    % colorbar;
    view(3)
end