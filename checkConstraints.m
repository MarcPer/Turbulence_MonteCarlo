clear; close all;
% Add dependencies
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
addpath('jsonlab');
addpath('propagation_routines');

% Get parameters
inputFile = ConfigHelper.getInputParametersFile;
simParams = SimulationParameters(loadjson(inputFile));

% Check constraints
simParams.constraintAnalysis;