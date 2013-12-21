%% Clear and close all
close all; clear all;
%
% This script simulates the effect of turbulence by modelling it as a
% series of phase screens. Technical details are taken from
%
% Numerical Simulation of Optical Wave Propagation
%            (with examples in MATLAB)
%               Jason D. Schmidt
%
% <LENGTH UNITS: m>

inOut = IOPaths;   % Imports functions and defines methods for export paths
usrIn = UserInput;

% Ask user to choose simulation type
usrIn.getSimulationType;

% Shutdown computer at the end of the script?
%   (0 = NO, 1 = SHUTDOWN, 2 = HIBERNATE)
usrIn.enumShutdown = 0;

%% SETUP AND SIMULATION PARAMETERS
% Setup geometry
simParams = SimulationParameters('FourthOrder', usrIn.isFourthOrder, 'Inverted', usrIn.isInverted);
constraintReturnCode = simParams.constraintAnalysis;

isAbort = usrIn.abortWhenConstraintFail(constraintReturnCode, usrIn.enumShutdown);
if isAbort
    return
end
close all; clear isAbort;

turbSimulator = TurbulenceSimulator(simParams);

%% RUN SIMULATION
intProfileGamma = turbSimulator.getIntensityForEachGamma('Normalized', true);

if turbSimulator.isAborted
    fprintf('Simulation aborted.');
    return;
end

%% COMPUTE RELEVANT VALUES
slitWidth = 50e-6;
[pwrSlit, plotInfo] = ...
    Calculator.computePowerThroughSlit(intProfileGamma, simParams);

%% PLOT RESULTS

% Beam transverse profiles at observation plane ('measured gammas only')
Plotter.plotIntensityProfilesForEachGamma(simParams.gammaStrength, intProfileGamma);

% Pslit = Pslit/Pslit(1);
% figure(2);
% set(gcf, 'Units', 'normalized', 'Position', [2/3 0.05 1/3 1/4]);
% plot(g, Pslit, 'LineWidth', 1.5);

% %% EXPORT RESULTS (only if simulation was completed)
% 
% if ~abort
%     % File name and folder
%     folder_save = fullfile(user_folder, 'Data', 'Computation', ...
%         'Phase_screen', date);
%     clock_str = clock;
%     date_str = [num2str(clock_str(1)), '-', num2str(clock_str(2)), ...
%         '-', num2str(clock_str(3))];
%     time_str = [num2str(clock_str(4)), ':', num2str(clock_str(5)), ...
%         ':', num2str(round(clock_str(6)))];
% 
%     if ~exist(folder_save, 'dir')
%         mkdir(folder_save)
%         curr_file = 1;
%     else
%         prev_files = dir(fullfile(folder_save, '*.dat'));
%         if isempty(prev_files)
%             curr_file = 1;
%         else
%             num_str = regexp(prev_files(end).name, '#(\d+).', 'tokens');
%             curr_file = str2num(num_str{1}{1});
%             curr_file = curr_file + 1;
%         end
%     end
% 
%     exp_file = [date_str, ' PSdata #', sprintf('%03u',curr_file),'.dat'];
%     fid = fopen(fullfile(folder_save,exp_file), 'w');
% 
%     % File info and model parameters
%     header = ['File saved at ', time_str, ' on the ', date_str];
%     fprintf(fid, '%s\n', header);
%     
%     fprintf(fid, 'Simulation for: %s\n', choice);
% 
%     data_params = char('Wavelength', 'Beam waist', 'Slit width', ...
%         'Inner scale', 'Outer scale', 'Grid spacing (source plane)', ...
%         'Grid spacing (obs. plane)', 'Grid size', ...
%         'Num. of realizations', 'Prop. distance', 'Num. of screens', ...
%         'Screen intervals', 'Distance between planes', ...
%         'Num. of prop. planes');
% 
%     fprintf(fid, '%s\n\n', '<< LENGTH UNITS: m >>');
%     format_spec = [data_params(1,:), '\t%3.2e\n', ...
%                 data_params(2,:), '\t%3.1e\n', ...
%                 data_params(3,:), '\t%3.1e\n', ...
%                 data_params(4,:), '\t%3.1e\n', ...
%                 data_params(5,:), '\t%3.1e\n', ...
%                 data_params(6,:), '\t%3.1e\n', ...
%                 data_params(7,:), '\t%3.1e\n', ...
%                 data_params(8,:), '\t%u\n', ...
%                 data_params(9,:), '\t%u\n', ...
%                 data_params(10,:), '\t%3.1e\n', ...
%                 data_params(11,:), '\t%u\n', ...
%                 data_params(12,:), '\t(%3.1e - %3.1e)\n', ...
%                 data_params(13,:), '\t%3.1e\n', ...
%                 data_params(14,:), '\t%u\n'];
%     
%     fprintf(fid, format_spec, wvl, wn, a0, ...
%         l0, L0, delta1, deltan, N, nreals, Dz, ...
%         nscr, zscr_min, zscr_max, mean(z(2:n)-z(1:n-1)),nz);
%     fprintf(fid, '\nPower over slit\n');
%     fprintf(fid, '%u\t', Pslit);
%     
%     % Close file
%     fclose(fid);
% end
% 
% %% SHUTDOWN COMPUTER?
% wtime = 90;         % Time limit to abort shutdown
% 
% if usrIn.enumShutdown
%     if usrIn.enumShutdown == 1
%         shut_mode = 'shutdown';
%         shut_cmd = 'shutdown /s';
%     else
%         shut_mode = 'hibernation';
%         shut_cmd = 'shutdown /h';
%     end
%     
%     fprintf('Computer %s enabled.\n', shut_mode);
%     tic;
%     reply_shut = 'k';
%     while (~sum(strcmpi(reply_shut, {'y', 'n', ''})) && toc < wtime)
%         if usrIn.enumShutdown == 1
%         reply_shut = waitinput('Proceed with shutdown? Y/N [Y]:    ', ...
%             wtime, 's');
%         else
%         reply_shut = waitinput('Proceed with hibernation? Y/N [Y]:    ', ...
%             wtime, 's');
%         end
%         fprintf('\n');
%     end
%     if (isempty(reply_shut) || ~sum(strcmpi(reply_shut, {'y', 'n'})))
%         reply_shut = 'Y';
%     end
%     if strcmpi(reply_shut,'y')
%         if usrIn.enumShutdown == 1
%             fprintf('\nShutting down.\n');
%         else
%             fprintf('\nHibernating.\n');
%         end
%         pause(3);
%         usrIn.enumShutdown = 0;
%         system(shut_cmd);
%     end
%     if usrIn.enumShutdown == 1
%         fprintf('Shutdown aborted.\n');
%     elseif usrIn.enumShutdown == 2
%         fprintf('Hibernation aborted.\n');
%     end
% end
