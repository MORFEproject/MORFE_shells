function SaveData(Cp, SSMParams, rdyn, NodeFile, ElementFile)
    str = NodeFile + ElementFile + "_" + ...
        num2str(SSMParams.max_order) + "_" + num2str(SSMParams.max_orderNA);
    folderName = "Output/" + str;
    MakeFolder(folderName);
    save(folderName+"\Cp.mat","Cp");
    save(folderName+"\rdyn.mat","rdyn");
end