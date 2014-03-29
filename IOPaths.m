classdef IOPaths<handle
    
    properties(GetAccess = public, SetAccess = private)
       dropboxFolder;
       inputParametersFileName;
       rootExportFolder;
       importFolders = {'Numerical simulation of optical wave propagation', ...
           'waitinput', 'hermite'};
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
        function expFileName = getExportFileName(io,data,varargin)
            ext = 'dat';
            identifier = 'SimData';
            if ~isempty(varargin)
                ext = varargin{1};
                if (length(varargin) > 1)
                    identifier = num2str(varargin{2});
                end
            end
            
            dateAndTime = regexp(data.info.date, '_', 'split');
            dt = dateAndTime{1};
            tm = dateAndTime{2};
            expFileName = [dt, '_', identifier, '_', tm, '.', ext];
            
            if exist(io.getExportPath, 'dir')
                files = io.getExportedFiles;
                fn = 0;
                if (ismember(expFileName,files))
                    fn = uigetfile(fullfile(io.getExportPath,['*.',ext]));
                end
                if fn
                    expFileName = fn;
                end
            end
        end
        function fData = openParametersFile(io)
            inputFolder = fullfile(io.dropboxFolder, 'MATLAB', 'Turbulence_MonteCarlo', 'inputParameters');
            fileName = fullfile(inputFolder, 'inputParameters.dat');
            fileName = uigetfile(fileName);
            fileName = fullfile(inputFolder, fileName);
            if (fileName == 0)
                error('IOPaths:noInputFile', 'No input parameters file selected');
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