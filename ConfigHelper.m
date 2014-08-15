classdef ConfigHelper
	methods(Static)
		function inFile = getInputParametersFile()
			disp('Choose input parameters file.');
			[inFolder, outFolder] = ConfigHelper.getIOFolders;
			[inFile, inFolder] = uigetfile(fullfile(inFolder,'*.json'), 'Choose input parameters file');
			if (inFile == 0)
				error('turbSimulator:noInputFile', 'No input file selected')
			end
			inFile = fullfile(inFolder, inFile);
			ConfigHelper.writeInputParametersFolder(inFolder, outFolder);
		end

		function setOutputFolder()
			disp('Select output folder for simulation results.');
			[inFolder, outFolder] = ConfigHelper.getIOFolders;
			if isempty(outFolder)
				outFolder = uigetdir('', 'Choose output folder for simulation results');
			end
			ConfigHelper.writeInputParametersFolder(inFolder, outFolder);
		end

		function [inFolder, outFolder] = getIOFolders()
			inFolder = '';
			outFolder = '';
			isNoConfigFile = isempty(ls('turbSimulatorConfig.json'));

			% Configuration file does not exist
			if (isNoConfigFile)
				return;
			end

			configStruct = loadjson('turbSimulatorConfig.json');
			if isfield(configStruct, 'inFolder')
				inFolder = configStruct.inFolder;
			end
			if isfield(configStruct, 'outFolder')
				outFolder = configStruct.outFolder;
			end
		end

		function writeInputParametersFolder(inFolder, outFolder)
			ioStruct = struct('inFolder', inFolder, 'outFolder', outFolder);
			savejson('', ioStruct, 'turbSimulatorConfig.json');
		end

	end

end
