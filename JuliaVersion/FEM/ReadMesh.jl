using DelimitedFiles
mutable struct PreCalculateParam
    K0::Array{Any, 1}
    M::Array{Any, 1}
    invKaa::Array{Any, 1}
    K0a::Array{Any, 1}
    Ka0::Array{Any, 1}
    det_a::Array{Float64, 2}
    
    D::Array{Any, 2}
    R::Array{Any, 2}
    G::Array{Any, 2}
    B0::Array{Any, 2}
    Bα::Array{Any, 2}
end

mutable struct MeshShell
    NodeCoord::Array{Float64, 2} 
    ElementNode::Array{Int64, 2}
    FIXDOF::Array{Int64, 1} 
    FREEDOF::Array{Int64, 1} 
    NEle::Int64 
    NNode::Int64 
    NDOF::Int64 
    SDOF::Array{Int64, 1} 
    NFIXDOF::Int64 
    NFREEDOF::Int64 
    ElementDOF::Array{Int64, 2} 
    ElementDOFNum::Int64 
    FREEDOF_list::Array{Int64, 1} 
    FixNode::Array{Int64, 1} 
    ElementNormal::Array{Any, 2} 
    PreCalculateParam::PreCalculateParam
    ForceNodeDOF::Array{Int64, 1}
    ForceVector::Array{Float64, 1}
    ForceNode::Array{Int64, 1}
end


function ReadShell(NodeFileName::String, ElementFileName::String)
    NodeData = readdlm(NodeFileName, skipstart=0)
    ElementData = readdlm(ElementFileName, skipstart=1)   
    NodeCoord = NodeData[:, 2:5]
    ElementNode = ElementData[:, 7:14] 
    ElementNode = convert(Array{Int64}, ElementNode) 
    MaxXDOF, MinXDOF, MaxYDOF, MinYDOF, MidNodeDOF, MaxXRow, MinXRow, MaxYRow, MinYRow, MidYRow, MidNodeRow, ForceNodeDOF, ForceNode = FindGeneralDofShell(NodeCoord)
    FixNode = [MaxXRow; MinXRow]
    FIXDOF = [MaxXDOF; MinXDOF]
    FIXDOF = vec(FIXDOF)
    FIXDOF = unique(FIXDOF)
    FIXDOF = sort(FIXDOF)
    NFIXDOF = length(FIXDOF)
    println("Fixed degree freedom = ", NFIXDOF, "\n")
    NNode = size(NodeCoord, 1)
    println("Total Node = ", NNode, "\n")
    NDOF = 6 * NNode
    println("Total degree freedom = ", NDOF, "\n")
    SDOF = collect(1:NDOF)
    FREEDOF = setdiff(SDOF, FIXDOF)
    NFREEDOF = length(FREEDOF)
    ForceVector = zeros(Float64, NFREEDOF)
    println("Free degree freedom = ", NFREEDOF, "\n")
    FREEDOF_list = 1:NFREEDOF
    NEle = size(ElementNode, 1)
    ElementDOF = zeros(Int64, (NEle, 48))
    for i = 1:NEle
        eleNode = ElementNode[i, :]
        for j = 1:8
            eleDOF = [eleNode[j] * 6 - 5, eleNode[j] * 6 - 4, eleNode[j] * 6 - 3, eleNode[j] * 6 - 2, eleNode[j] * 6 - 1, eleNode[j] * 6]
            ElementDOF[i, 6 * j - 5:6 * j] = eleDOF
        end
    end
    ElementDOFNum = 48
    # =============== calculate element normal ===============
    ElementNormal = Array{Any}(undef, NEle, 8)
    for i = 1 :  NEle
        eleNode = ElementNode[i, :]
        x1 = NodeCoord[eleNode[1], 1]
        y1 = NodeCoord[eleNode[1], 2]
        z1 = NodeCoord[eleNode[1], 3]
        x2 = NodeCoord[eleNode[2], 1]
        y2 = NodeCoord[eleNode[2], 2]
        z2 = NodeCoord[eleNode[2], 3]
        x3 = NodeCoord[eleNode[3], 1]
        y3 = NodeCoord[eleNode[3], 2]
        z3 = NodeCoord[eleNode[3], 3]
        x4 = NodeCoord[eleNode[4], 1]
        y4 = NodeCoord[eleNode[4], 2]
        z4 = NodeCoord[eleNode[4], 3]
        x5 = NodeCoord[eleNode[5], 1]
        y5 = NodeCoord[eleNode[5], 2]
        z5 = NodeCoord[eleNode[5], 3]
        x6 = NodeCoord[eleNode[6], 1]
        y6 = NodeCoord[eleNode[6], 2]
        z6 = NodeCoord[eleNode[6], 3]
        x7 = NodeCoord[eleNode[7], 1]
        y7 = NodeCoord[eleNode[7], 2]
        z7 = NodeCoord[eleNode[7], 3]
        x8 = NodeCoord[eleNode[8], 1]
        y8 = NodeCoord[eleNode[8], 2]
        z8 = NodeCoord[eleNode[8], 3]

        xmid = (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8) / 8
        ymid = (y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8) / 8
        zmid = (z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) / 8

        vector1 = [x1, y1, z1]
        vector2 = [x2, y2, z2]
        vector3 = [x3, y3, z3]
        vector4 = [x4, y4, z4]
        vector5 = [x5, y5, z5]
        vector6 = [x6, y6, z6]
        vector7 = [x7, y7, z7]
        vector8 = [x8, y8, z8]
        vectormid = [xmid, ymid, zmid]

        e1 = cross(vector5-vector1, vector8-vector1)
        e2 = cross(vector6-vector2, vector5-vector2)
        e3 = cross(vector7-vector3, vector6-vector3)
        e4 = cross(vector8-vector4, vector7-vector4)

        e5 = cross(vector2-vector5, vectormid-vector5)
        e6 = cross(vector3-vector6, vectormid-vector6)
        e7 = cross(vector4-vector7, vectormid-vector7)
        e8 = cross(vector1-vector8, vectormid-vector8)

        ElementNormal[i, 1] = e1/norm(e1)
        ElementNormal[i, 2] = e2/norm(e2)
        ElementNormal[i, 3] = e3/norm(e3)
        ElementNormal[i, 4] = e4/norm(e4)
        ElementNormal[i, 5] = e5/norm(e5)
        ElementNormal[i, 6] = e6/norm(e6)
        ElementNormal[i, 7] = e7/norm(e7)
        ElementNormal[i, 8] = e8/norm(e8)

    end
    # average normal calculation
    ShareNodeCell = Array{Any}(undef, NNode, 2)
    for NodeNum = 1 : NNode
        indices = findall(x -> x == NodeNum, ElementNode)
        rows = [index[1] for index in indices]
        cols = [index[2] for index in indices]
        ShareNodeCell[NodeNum, 1] = [rows cols]
        ShareNodeCell[NodeNum, 2] = length(rows)
    end
    NodeNormalDir = zeros(Float64, (NNode, 3))
    for NodeNum = 1 : NNode
        indices = ShareNodeCell[NodeNum, 1]
        for i = 1 : ShareNodeCell[NodeNum, 2]
            ElementNum = indices[i, 1]
            ElementNodeNum = indices[i, 2]
            NodeNormalDir[NodeNum, :] += ElementNormal[ElementNum, ElementNodeNum]
        end
        NodeNormalDir[NodeNum, :] = NodeNormalDir[NodeNum, :] / ShareNodeCell[NodeNum, 2]
    end
    for i = 1 : NEle
        eleNode = ElementNode[i, :]
        for j = 1 : 8
            ElementNormal[i, j] = NodeNormalDir[eleNode[j], :]
        end
    end
    # ============ end of average normal calculation ==============
    # initialization 
    K0 = Array{Any}(undef, NEle)
    M = Array{Any}(undef, NEle)
    invKaa = Array{Any}(undef, NEle)
    K0a = Array{Any}(undef, NEle)
    Ka0 = Array{Any}(undef, NEle)
    det_a = zeros(Int64, (NEle, 1))
    D = Array{Any}(undef, NEle, 8)
    R = Array{Any}(undef, NEle, 8)
    G = Array{Any}(undef, NEle, 8)
    B0 = Array{Any}(undef, NEle, 8)
    Bα = Array{Any}(undef, NEle, 8)
    PreParam = PreCalculateParam(K0, M, invKaa, K0a, Ka0, det_a, D, R, G, B0, Bα)
    return MeshShell(NodeCoord, ElementNode, FIXDOF, FREEDOF, NEle, NNode, NDOF, SDOF, NFIXDOF, NFREEDOF, ElementDOF, ElementDOFNum, FREEDOF_list, FixNode, ElementNormal,PreParam,ForceNodeDOF,ForceVector,ForceNode)

end


function addThicknessColumn(inputFile::String, outputFile::String, tmin, tmax)
    data = readdlm(inputFile)
    nodeID = data[:, 1]
    xCoord = data[:, 2]
    yCoord = data[:, 3]
    zCoord = data[:, 4]
    
    xMin = minimum(xCoord)
    xMax = maximum(xCoord)
    
    #   t(x) = tmin + (tmax - tmin)*((x - xMin)/(xMax - xMin))
    thickness = tmin .+ (tmax - tmin) .* ((xCoord .- xMin) ./ (xMax - xMin))
    
    outData = hcat(nodeID, xCoord, yCoord, zCoord, thickness)
    writedlm(outputFile, outData, ' ')
    println("New node document with thickness: $outputFile")
end


function FindAllFixedDofShell(Node::Vector{Int64})
    DofList = [Node .* 6 .- 5; Node .* 6 .- 4; Node .* 6 .- 3; Node .* 6 .- 2; Node .* 6 .- 1; Node .* 6]
    DofList = reshape(DofList, (6 * length(Node), 1))
    return DofList    
end

function FindSSSSDofShell(Node::Vector{Int64})
    DofList = [Node .* 6 .- 5; Node .* 6 .- 4; Node .* 6 .- 3;]
    DofList = reshape(DofList, (3 * length(Node), 1))
    return DofList    
end


function FindGeneralDofShell(NodeCoord::Array{Float64, 2})
    MaxX = maximum(NodeCoord[:, 1])
    MinX = minimum(NodeCoord[:, 1])
    MaxY = maximum(NodeCoord[:, 2])
    MinY = minimum(NodeCoord[:, 2])
    MaxZ = maximum(NodeCoord[:, 3])
    MinZ = minimum(NodeCoord[:, 3])
    MaxXRow = findall(x -> x == MaxX, NodeCoord[:, 1])
    MinXRow = findall(x -> x == MinX, NodeCoord[:, 1])
    MidX = (MaxX + MinX) / 2
    MidXRow = findall(x -> x == MidX, NodeCoord[:, 1])
    MaxYRow = findall(x -> x == MaxY, NodeCoord[:, 2])
    MinYRow = findall(x -> x == MinY, NodeCoord[:, 2])
    MidY = (MaxY + MinY) / 2
    MidYRow = findall(x -> x == MidY, NodeCoord[:, 2])
    MidNodeRow = intersect(MidYRow, MidXRow)
    MaxZRow = findall(x -> x == MaxZ, NodeCoord[:, 3])
    MinZRow = findall(x -> x == MinZ, NodeCoord[:, 3])
    MaxXDOF = FindSSSSDofShell(MaxXRow)
    MinXDOF = FindSSSSDofShell(MinXRow)
    MaxYDOF = FindSSSSDofShell(MaxYRow)
    MinYDOF = FindSSSSDofShell(MinYRow)
    MidNodeDOF = FindAllFixedDofShell(MidNodeRow)
    ForceNode = findall((abs.(NodeCoord[:, 1] .+ (MaxX - MinX) / 4) .<= 1e-3) .& (NodeCoord[:, 2] .== (MaxY + MinY) / 2))
    ForceNodeDOF = ForceNode .* 6 .- 3
    return MaxXDOF, MinXDOF, MaxYDOF, MinYDOF, MidNodeDOF, MaxXRow, MinXRow, MaxYRow, MinYRow, MidYRow, MidNodeRow, ForceNodeDOF, ForceNode

end