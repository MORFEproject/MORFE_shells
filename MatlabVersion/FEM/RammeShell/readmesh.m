function [Mesh] = readmesh(NodeFile,ElementFile, BoundaryType, DirectionType)
    node = textread(NodeFile);
    element = textread(ElementFile);
    element = element(:,[7:14]);
    node = node(:,[2:5]);
    xmin = min(node(:,1));
    xmax = max(node(:,1));
    xmid = (xmin+xmax)/2;

    ymin = min(node(:,2));
    ymax = max(node(:,2));
    ymid = (ymin+ymax)/2;

    zmin = min(node(:,3));
    zmax = max(node(:,3));
    zmid = (zmin+zmax)/2;    
    fixnode1 = find(node(:,1)==xmin);
    fixnode2 = find(node(:,1)==xmax);
    fixnode3 = find(node(:,2)==ymin);
    fixnode4 = find(node(:,2)==ymax);   
    fixnodez = find(node(:,3)==zmax);

    fixnode_core1 = find(node(:,1)==xmin & node(:,2)==ymin);
    fixnode_core2 = find(node(:,1)==xmin & node(:,2)==ymax);
    fixnode_core3 = find(node(:,1)==xmax & node(:,2)==ymin);
    fixnode_core4 = find(node(:,1)==xmax & node(:,2)==ymax);

    forcenode = find(abs(node(:,1)+(xmax-xmin)/4)<=1e-3 & node(:,2)==ymid);
    forcedof = forcenode*6-3;
    switch BoundaryType
        case 'c' % 
            fixnode = [fixnode1];
            fixdof1 = [6*fixnode1-5, 6*fixnode1-4, 6*fixnode1-3, 6*fixnode1-2, 6*fixnode1-1, 6*fixnode1];
            fixdof = [fixdof1(:)];
            fixdof = unique(fixdof);
            fixdof = sort(fixdof);            
        case 'c-c' % 
            fixnode = [fixnode1; fixnode2];
            fixdof1 = [6*fixnode1-5, 6*fixnode1-4, 6*fixnode1-3, 6*fixnode1-2, 6*fixnode1-1, 6*fixnode1];
            fixdof2 = [6*fixnode2-5, 6*fixnode2-4, 6*fixnode2-3, 6*fixnode2-2, 6*fixnode2-1, 6*fixnode2];
            fixdof = [fixdof1(:); fixdof2(:)];
            fixdof = unique(fixdof);
            fixdof = sort(fixdof);
        case 's-s' % 
            fixnode = [fixnode1; fixnode2];
            fixdof1 = [6*fixnode1-5, 6*fixnode1-4, 6*fixnode1-3];
            fixdof2 = [6*fixnode2-5, 6*fixnode2-4, 6*fixnode2-3];
            fixdof = [fixdof1, fixdof2];
            fixdof = unique(fixdof);
            fixdof = sort(fixdof);            
        case 'c-c-c-c'
            fixnode = [fixnode1; fixnode2; fixnode3; fixnode4];
            fixdof1 = [6*fixnode1-5, 6*fixnode1-4, 6*fixnode1-3, 6*fixnode1-2, 6*fixnode1-1, 6*fixnode1];
            fixdof2 = [6*fixnode2-5, 6*fixnode2-4, 6*fixnode2-3, 6*fixnode2-2, 6*fixnode2-1, 6*fixnode2];
            fixdof3 = [6*fixnode3-5, 6*fixnode3-4, 6*fixnode3-3, 6*fixnode3-2, 6*fixnode3-1, 6*fixnode3];
            fixdof4 = [6*fixnode4-5, 6*fixnode4-4, 6*fixnode4-3, 6*fixnode4-2, 6*fixnode4-1, 6*fixnode4];
            fixdof = [fixdof1(:); fixdof2(:); fixdof3(:); fixdof4(:)];
            fixdof = unique(fixdof);
            fixdof = sort(fixdof);
        case 's-s-s-s'
            fixnode = [fixnode1; fixnode2; fixnode3; fixnode4];
            fixdof1 = [6*fixnode1-5, 6*fixnode1-4, 6*fixnode1-3];
            fixdof2 = [6*fixnode2-5, 6*fixnode2-4, 6*fixnode2-3];
            fixdof3 = [6*fixnode3-5, 6*fixnode3-4, 6*fixnode3-3];
            fixdof4 = [6*fixnode4-5, 6*fixnode4-4, 6*fixnode4-3];
            fixdof = [fixdof1(:); fixdof2(:); fixdof3(:); fixdof4(:)];
            fixdof = unique(fixdof);
            fixdof = sort(fixdof);
    end
    % 
    ElementDof = zeros(size(element,1),48);
    for i = 1 : size(element,1)
        for j = 1 : 8
            ElementDof(i,6*j-5) = 6*element(i,j)-5;
            ElementDof(i,6*j-4) = 6*element(i,j)-4;
            ElementDof(i,6*j-3) = 6*element(i,j)-3;
            ElementDof(i,6*j-2) = 6*element(i,j)-2;
            ElementDof(i,6*j-1) = 6*element(i,j)-1;
            ElementDof(i,6*j) = 6*element(i,j);
        end
    end
    switch DirectionType
        case 'Independent'
            % solve normal of element
            ElementDirection = cell(size(element,1),8);
            for i = 1 : size(element,1)
                n1 = element(i,1); x1 = node(n1,1);  y1 = node(n1,2); z1 = node(n1,3);
                n2 = element(i,2); x2 = node(n2,1);  y2 = node(n2,2); z2 = node(n2,3);
                n3 = element(i,3); x3 = node(n3,1);  y3 = node(n3,2); z3 = node(n3,3);
                n4 = element(i,4); x4 = node(n4,1);  y4 = node(n4,2); z4 = node(n4,3);
                n5 = element(i,5); x5 = node(n5,1);  y5 = node(n5,2); z5 = node(n5,3);
                n6 = element(i,6); x6 = node(n6,1);  y6 = node(n6,2); z6 = node(n6,3);
                n7 = element(i,7); x7 = node(n7,1);  y7 = node(n7,2); z7 = node(n7,3);
                n8 = element(i,8); x8 = node(n8,1);  y8 = node(n8,2); z8 = node(n8,3);
                xmid = (x1+x2+x3+x4+x5+x6+x7+x8)/8;
                ymid = (y1+y2+y3+y4+y5+y6+y7+y8)/8;
                zmid = (z1+z2+z3+z4+z5+z6+z7+z8)/8;
        
                vector1 = [x1, y1, z1]; vector2 = [x2, y2, z2]; vector3 = [x3, y3, z3]; 
                vector4 = [x4, y4, z4]; vector5 = [x5, y5, z5]; vector6 = [x6, y6, z6];
                vector7 = [x7, y7, z7]; vector8 = [x8, y8, z8]; vectormid = [xmid, ymid, zmid];
        
                e1 = cross((vector5-vector1),(vector8-vector1));
                e2 = cross((vector6-vector2),(vector5-vector2));
                e3 = cross((vector7-vector3),(vector6-vector3));
                e4 = cross((vector8-vector4),(vector7-vector4));
        
                e5 = cross((vector2-vector5),(vectormid-vector5));
                e6 = cross((vector3-vector6),(vectormid-vector6));
                e7 = cross((vector4-vector7),(vectormid-vector7));
                e8 = cross((vector1-vector8),(vectormid-vector8));
                
                ElementDirection{i,1} = e1/norm(e1);
                ElementDirection{i,2} = e2/norm(e2);
                ElementDirection{i,3} = e3/norm(e3);
                ElementDirection{i,4} = e4/norm(e4);
                ElementDirection{i,5} = e5/norm(e5);
                ElementDirection{i,6} = e6/norm(e6);
                ElementDirection{i,7} = e7/norm(e7);
                ElementDirection{i,8} = e8/norm(e8);
        
            end    
        case 'Average'
            % solve normal of element
            ElementDirection = cell(size(element,1),8);
            for i = 1 : size(element,1)
                n1 = element(i,1); x1 = node(n1,1);  y1 = node(n1,2); z1 = node(n1,3);
                n2 = element(i,2); x2 = node(n2,1);  y2 = node(n2,2); z2 = node(n2,3);
                n3 = element(i,3); x3 = node(n3,1);  y3 = node(n3,2); z3 = node(n3,3);
                n4 = element(i,4); x4 = node(n4,1);  y4 = node(n4,2); z4 = node(n4,3);
                n5 = element(i,5); x5 = node(n5,1);  y5 = node(n5,2); z5 = node(n5,3);
                n6 = element(i,6); x6 = node(n6,1);  y6 = node(n6,2); z6 = node(n6,3);
                n7 = element(i,7); x7 = node(n7,1);  y7 = node(n7,2); z7 = node(n7,3);
                n8 = element(i,8); x8 = node(n8,1);  y8 = node(n8,2); z8 = node(n8,3);
                xmid = (x1+x2+x3+x4+x5+x6+x7+x8)/8;
                ymid = (y1+y2+y3+y4+y5+y6+y7+y8)/8;
                zmid = (z1+z2+z3+z4+z5+z6+z7+z8)/8;
        
                vector1 = [x1, y1, z1]; vector2 = [x2, y2, z2]; vector3 = [x3, y3, z3]; 
                vector4 = [x4, y4, z4]; vector5 = [x5, y5, z5]; vector6 = [x6, y6, z6];
                vector7 = [x7, y7, z7]; vector8 = [x8, y8, z8]; vectormid = [xmid, ymid, zmid];
        
                e1 = cross((vector5-vector1),(vector8-vector1));
                e2 = cross((vector6-vector2),(vector5-vector2));
                e3 = cross((vector7-vector3),(vector6-vector3));
                e4 = cross((vector8-vector4),(vector7-vector4));
        
                e5 = cross((vector2-vector5),(vectormid-vector5));
                e6 = cross((vector3-vector6),(vectormid-vector6));
                e7 = cross((vector4-vector7),(vectormid-vector7));
                e8 = cross((vector1-vector8),(vectormid-vector8));
                
                ElementDirection{i,1} = e1/norm(e1);
                ElementDirection{i,2} = e2/norm(e2);
                ElementDirection{i,3} = e3/norm(e3);
                ElementDirection{i,4} = e4/norm(e4);
                ElementDirection{i,5} = e5/norm(e5);
                ElementDirection{i,6} = e6/norm(e6);
                ElementDirection{i,7} = e7/norm(e7);
                ElementDirection{i,8} = e8/norm(e8);
        
            end 
            % =============== average norm =================
            ShareNodeCell = cell(size(node, 1), 2);
            for NodeNum = 1 : size(node, 1)
                [row, col] = find(element == NodeNum);
                ShareNodeCell{NodeNum, 1} = [row, col];
                ShareNodeCell{NodeNum, 2} = length(row);
            end
            % 
            NodeNormalDir = zeros(size(node, 1), 3);
            for i = 1 : size(node, 1)
                NormalDir = zeros(1, 3);
                ShareNodeMatrix = ShareNodeCell{i, 1};
                for j = 1 : ShareNodeCell{i, 2} 
                    NormalDir = NormalDir +  ElementDirection{ShareNodeMatrix(j,1), ShareNodeMatrix(j,2)};         
                end
                NodeNormalDir(i, :) = NormalDir / ShareNodeCell{i, 2};
            end
            % 
            for i = 1 : size(element,1)
                n1 = element(i,1); 
                n2 = element(i,2);
                n3 = element(i,3);
                n4 = element(i,4);
                n5 = element(i,5);
                n6 = element(i,6);
                n7 = element(i,7);
                n8 = element(i,8);
        
                ElementDirection{i,1} = NodeNormalDir(n1,:); ElementDirection{i,2} = NodeNormalDir(n2,:);
                ElementDirection{i,3} = NodeNormalDir(n3,:); ElementDirection{i,4} = NodeNormalDir(n4,:);
                ElementDirection{i,5} = NodeNormalDir(n5,:); ElementDirection{i,6} = NodeNormalDir(n6,:);
                ElementDirection{i,7} = NodeNormalDir(n7,:); ElementDirection{i,8} = NodeNormalDir(n8,:);
        
            end
            % =============== end average norm =================
    end

    sdof = 6*size(node,1);
    sdoflist = 1:sdof;
    freedof = setdiff(sdoflist, fixdof);
    TotalNode = 1 : size(node,1);
    freenode = setdiff(TotalNode, fixnode);

    % data package
    Mesh.Node = node;
    Mesh.NodeNum = size(node,1);
    Mesh.ElementDof = ElementDof;
    Mesh.ElementDirection = ElementDirection;
    Mesh.Element = element;
    Mesh.Freedof = freedof;
    Mesh.FreedofNum = length(freedof);
    Mesh.Sdof = sdof;
    Mesh.ElementNum = size(element,1);
    Mesh.NodeNum = size(node,1);
    Mesh.Fixdof = fixdof;
    Mesh.FixdofNum = length(fixdof);
    Mesh.FixNode = fixnode;
    Mesh.FreeNode = freenode;
    Mesh.ForceNode = forcenode;
    Mesh.ForceDof = forcedof;
end