function PlotMesh(Mesh)
    element = Mesh.Element;
    node = Mesh.Node;
    fixnode = Mesh.FixNode;
    forcenode = Mesh.ForceNode;
    freenode = Mesh.FreeNode;
    elementum = Mesh.ElementNum;
    face = [1 2 3 4;];
    figure
    hold on 
    for i = 1:elementum
        vert=[node(element(i,1),1) node(element(i,1),2) node(element(i,1),3);
              node(element(i,2),1) node(element(i,2),2) node(element(i,2),3);
              node(element(i,3),1) node(element(i,3),2) node(element(i,3),3);
              node(element(i,4),1) node(element(i,4),2) node(element(i,4),3)];
        patch('Faces',face,'Vertices',vert,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',1);  
    end
    for i = 1 : length(fixnode)
        x = node(fixnode(i),1); y = node(fixnode(i),2); z = node(fixnode(i),3);
        scatter3(x,y,z,'filled','red')
    end    
    for i = 1 : length(forcenode)
        x = node(forcenode(i),1); y = node(forcenode(i),2); z = node(forcenode(i),3);
        scatter3(x,y,z,'filled','yellow')
    end  
    axis equal; 
    axis off; 
    view(3)
end