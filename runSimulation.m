clear; close all;

% Add dependencies
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
addpath('jsonlab');
addpath('propagation_routines');
addpath('waitinput');

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
	ticID = tic;
	results = turbSimulator.simulate;
catch exception
    UserInput.shutdownComputer(shutdownOrHibernate);
    set(0,'ShowHiddenHandles','on')
    delete(get(0,'Children'))
    rethrow(exception);
end

% Print elapsed time in readable form
tm = toc(ticID);
disp(Util.printReadableTime(tm));

% Export raw results
fileName = Exporter.exportToDisk(results, simParams, inputFile);

% Plot results
plotter = Plotter(results.params.simulationType);
plotter.plot(results);

% Send email with results (if configured)
sendEmail({[fileName '.dat'], [fileName '.mat']});

% Save irradiance plots to disk (disabled as figure might take considerable disk space)
% Exporter.exportFigure(results, simParams.simulationType)

% Shutdown computer
UserInput.shutdownComputer(shutdownOrHibernate);