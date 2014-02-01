%% Clear and close all
close all; clear all;
%
% This script simulates the effect of turbulence, modelling it as a
% series of phase screens. Technical details are taken from
%
% Numerical Simulation of Optical Wave Propagation
%            (with examples in MATLAB)
%               Jason D. Schmidt
%
% <LENGTH UNITS: m>

ioPaths = IOPaths;   % Imports functions and defines methods for export paths
usrIn = UserInput;

% Ask user to choose simulation type
usrIn.getSimulationType;

% Shutdown computer at the end of the script?
%   (0 = NO, 1 = SHUTDOWN, 2 = HIBERNATE)
usrIn.enumShutdown = 2;

%% SETUP AND SIMULATION PARAMETERS
% Setup geometry
simParams = SimulationParameters(ioPaths,'FourthOrder', usrIn.isFourthOrder, 'Inverted', usrIn.isInverted);
simParams.setPointDetectorAtObservationPlane(true);
constraintReturnCode = simParams.constraintAnalysis;

isAbort = usrIn.abortWhenConstraintFail(constraintReturnCode, usrIn.enumShutdown);
if isAbort
    return
end
close all; clear isAbort;

turbSimulator = TurbulenceSimulator(simParams);

%% RUN SIMULATION
try
    apertureRadius = simParams.waistAtObservationPlane/2;
    pwrGamma = turbSimulator.getPowerOnCircularApertureForEachGamma(apertureRadius,'Normalized', true);
    %pwrGamma = turbSimulator.getIrradianceForEachGamma('Normalized', true);
catch exception
    % %% SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
    rethrow(exception);
end

if turbSimulator.isAborted
    fprintf('Simulation aborted.\n');
    return;
end

%% PLOT RESULTS

% Beam transverse profiles at observation plane ('measured gammas only')
close all;

try
    %Plotter.plotIntensityProfilesForEachGamma(pwrGamma);
    Plotter.plot2D(pwrGamma)
    
    % EXPORT RESULTS (only if simulation was completed)
    Exporter.exportToDisk(ioPaths, pwrGamma, usrIn, turbSimulator.simulationParameters);
catch exception
    % SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
    rethrow(exception);
end
usrIn.shutdownComputer;