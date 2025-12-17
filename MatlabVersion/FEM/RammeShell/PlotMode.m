function PlotMode(Mesh, MatParams, SSMParams, ModeNum, Scale)
    element = Mesh.Element;
    node = Mesh.Node;
    elementnum = Mesh.ElementNum;
    face = [1 2 3 4];
    [K,M] = AssembleKM(Mesh, MatParams);
    [Phi, Lambda] = eigs(K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),SSMParams.ComputeMode,'sm');
    Phi1 = Phi(:,ModeNum);
    if sum(Phi1)<0
       Phi1 = -Phi1;
    end
    PhiFull = zeros(Mesh.Sdof,1);
    PhiFull(Mesh.Freedof) = Phi1;
    PhiU = reshape(PhiFull,6,[])';
    figure
    title(['Mode' num2str(ModeNum)])
    hold on 
    for i = 1 : elementnum
        Ue = PhiU(element(i,:),:);
        Ue_plane = Ue(:, [1:3]);
        we = Ue(:,3);
        wreduce = we(face);
        NodeLocal = node(element(i,:),:) + Scale * Ue_plane;
        NodeLocal2 = node(element(i,:),:);
        vert = NodeLocal(face,:);
        vert2 = NodeLocal2(face,:);
        patch(vert2(:,1),vert2(:,2),vert2(:,3),[0 0.4470 0.7410],'FaceAlpha',0.2);
        patch(vert(:,1),vert(:,2),vert(:,3),wreduce,'FaceColor','interp');
    end
    axis equal; axis off; box off;
    % colorbar;
    view(3)
end