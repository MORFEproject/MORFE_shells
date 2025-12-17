function InitModule(varargin)
    warning('off')
    numOptional = numel(varargin);
    thisFolder = fileparts(mfilename('fullpath'));
    removeLocalPaths(thisFolder);
    
    if numOptional == 0
        targetFolderFullPath = fullfile(thisFolder, 'Tool');
        addpath(genpath(targetFolderFullPath));    
    
        targetFolderFullPath = fullfile(thisFolder, 'Output');
        addpath(genpath(targetFolderFullPath));    
    else
        targetFolderFullPath = fullfile(thisFolder, 'Tool');
        addpath(genpath(targetFolderFullPath));    
    
        targetFolderFullPath = fullfile(thisFolder, 'Output');
        addpath(genpath(targetFolderFullPath));    

        FEMType = varargin{1};
        AnalysisType = varargin{2};
        switch FEMType
            case 'RammeShell'
                addModulePath(thisFolder, 'FEM/RammeShell');
                addModulePath(thisFolder, 'Material/RammeShell');
            case 'DPIM'
                addModulePath(thisFolder, 'FEM/DPIM');
                addModulePath(thisFolder, 'Material/DPIM');
        end
    
        switch AnalysisType
            case 'DPIM'
                addModulePath(thisFolder, 'DPIM');
                addModulePath(thisFolder, 'FEM/FEMDPIM');
            case 'HB'
                addModulePath(thisFolder, 'HB');
        end

    end

    function addModulePath(rootFolder, subFolder)
        targetFolder = fullfile(rootFolder, subFolder);
        if isfolder(targetFolder)
            addpath(genpath(targetFolder));
        else
            warning('The folder %s does not exist.', targetFolder);
        end
    end

    function removeLocalPaths(rootFolder)
        localPaths = genpath(rootFolder);
        folderList = strsplit(localPaths, pathsep);
        for k = 1:length(folderList)
            if ~isempty(folderList{k}) && isfolder(folderList{k})
                rmpath(folderList{k});
            end
        end
    end
end
