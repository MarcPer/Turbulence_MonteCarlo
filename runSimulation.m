clear; close all;

% Add dependencies
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
addpath('jsonlab');
addpath('propagation_routines');

% Setup simulation
inputFile = ConfigHelper.getInputParametersFile;
ConfigHelper.setOutputFolder;
simParams = SimulationParameters(loadjson(inputFile));
shutdownOrHibernate = UserInput.confirmShutdownOrHibernate(simParams.shutdownOrHibernate);

% Check constraints
isAbort = simParams.checkConstraints;
if isAbort
	return;
end

% Perform simulation
turbSimulator = TurbulenceSimulator(simParams);
try
	results = turbSimulator.simulate;
catch exception
    UserInput.shutdownComputer(shutdownOrHibernate);
    set(0,'ShowHiddenHandles','on')
    delete(get(0,'Children'))
    rethrow(exception);
end

% Export results
Exporter.exportToDisk(results, simParams, inputFile);


% Plot results
plotter = Plotter(simParams.simulationType);
plotter.plot(results);

% Shutdown computer
UserInput.shutdownComputer(shutdownOrHibernate);