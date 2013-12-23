classdef IOPaths<handle
    
    properties(GetAccess = public, SetAccess = private)
       dropboxFolder;
       inputParametersFileName;
       rootExportFolder;
       importFolders = {'Numerical simulation of optical wave propagation', ...
           'waitinput'};
    end
    
    methods (Access = public)
        function io = IOPaths()
            usrFolder = userpath;
            usrFolder(end) = [];
            io.dropboxFolder = fullfile(usrFolder, '../..', 'Dropbox');
            
            % Import folders
            for impt = 1 : numel(io.importFolders)
                addpath(fullfile('./', io.importFolders{impt}));
            end
            
            io.rootExportFolder = fullfile(io.dropboxFolder, 'Data', ...
                'Computation', 'Phase_screen');
        end
        function expPath = getExportPath(io)
           dt =  datestr(date, 'yyyy-mm-dd');
           expPath = fullfile(io.rootExportFolder,dt);
        end
        function expFileName = getExportFileName(io)
            dt = datestr(date, 'yyyy-mm-dd');
            tm = datestr(clock, 'HHMMSS');
            expFileName = [dt, '_SimData_', tm, '.dat'];
            
            if exist(io.getExportPath, 'dir')
                files = io.getExportedFiles;
                fn = 0;
                if (ismember(expFileName,files))
                    fn = uigetfile(fullfile(io.getExportPath,'*.dat'));
                end
                if fn
                    expFileName = fn;
                end
            end
        end
        function fData = openParametersFile(io)
            fileName = 'inputParameters.dat';
            fileName = uigetfile(fileName);
            if (fileName == 0)
                return;
            end
            io.inputParametersFileName = fileName;
            fid = fopen(fileName);
            try
                fData = fread(fid, inf, '*char');
                fData = fData';
            catch exception
                fprintf('Error reading %s.\n',fileName);
                fclose(fid);
                rethrow(exception);
            end
            fclose(fid);
        end
    end
    
    methods(Access = private)
        function fileNames = getExportedFiles(io)
            fileNamesChar = ls(io.getExportPath);
            [numFiles, ~] = size(fileNamesChar);
            fileNames = cell(numFiles,1);
            
            for f = 1 : numFiles
                fileNames{f} = strtrim(fileNamesChar(f,:));
            end
        end
    end
end