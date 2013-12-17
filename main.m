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
Phase_screen_choose_simulation;

% Shutdown computer at the end of the script?
%   (0 = NO, 1 = SHUTDOWN, 2 = HIBERNATE)
shut = 0;
Phase_screen_shutdown_input;

%% SETUP AND SIMULATION PARAMETERS
% Setup geometry
simParams = SimulationParameters('FourthOrder', true);
constraintReturnCode = simParams.constraintAnalysis;

isAbort = UserInput.abortWhenConstraintFail(constraintReturnCode, shut);
if isAbort
    return
end

turbSimulator = TurbulenceSimulator(simParams);

%% RUN SIMULATION
intensityProfileForEachTurbStrength = ...
    turbSimulator.getIntensityForEachGamma();

%% PLOT RESULTS (if computation was not aborted)

if ~abort
% Beam transverse profiles at observation plane ('measured gammas only')
figure(1);
set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);

for g0i = 1 : leng0
    subplot( ceil(sqrt(leng0)), ceil(sqrt(leng0)), g0i);
    if ~isempty(find(Iout(:, :, g0indx(g0i)),1))
        imagesc(Iout(:, :, g0indx(g0i)));
        axis off;
        title(['\gamma = ', num2str(g0(g0i)*1e6), ' \mum']);
    end
end

% Power going through slit
if coinc
    for idxStrength = 1 : leng
    A = Iout(:,:,idxStrength);
    Pslit(idxStrength) = CoincSlitIntegrate(A,a);
    end
else
    for idxStrength = 1 : leng
        A = Iout(:,:,idxStrength);
        Pslit(idxStrength) = SlitIntegrate(A,a);
    end
end

Pslit = Pslit/Pslit(1);
figure(2);
set(gcf, 'Units', 'normalized', 'Position', [2/3 0.05 1/3 1/4]);
plot(g, Pslit, 'LineWidth', 1.5);
end

%% EXPORT RESULTS (only if simulation was completed)

if ~abort
    % File name and folder
    folder_save = fullfile(user_folder, 'Data', 'Computation', ...
        'Phase_screen', date);
    clock_str = clock;
    date_str = [num2str(clock_str(1)), '-', num2str(clock_str(2)), ...
        '-', num2str(clock_str(3))];
    time_str = [num2str(clock_str(4)), ':', num2str(clock_str(5)), ...
        ':', num2str(round(clock_str(6)))];

    if ~exist(folder_save, 'dir')
        mkdir(folder_save)
        curr_file = 1;
    else
        prev_files = dir(fullfile(folder_save, '*.dat'));
        if isempty(prev_files)
            curr_file = 1;
        else
            num_str = regexp(prev_files(end).name, '#(\d+).', 'tokens');
            curr_file = str2num(num_str{1}{1});
            curr_file = curr_file + 1;
        end
    end

    exp_file = [date_str, ' PSdata #', sprintf('%03u',curr_file),'.dat'];
    fid = fopen(fullfile(folder_save,exp_file), 'w');

    % File info and model parameters
    header = ['File saved at ', time_str, ' on the ', date_str];
    fprintf(fid, '%s\n', header);
    
    fprintf(fid, 'Simulation for: %s\n', choice);

    data_params = char('Wavelength', 'Beam waist', 'Slit width', ...
        'Inner scale', 'Outer scale', 'Grid spacing (source plane)', ...
        'Grid spacing (obs. plane)', 'Grid size', ...
        'Num. of realizations', 'Prop. distance', 'Num. of screens', ...
        'Screen intervals', 'Distance between planes', ...
        'Num. of prop. planes');

    fprintf(fid, '%s\n\n', '<< LENGTH UNITS: m >>');
    format_spec = [data_params(1,:), '\t%3.2e\n', ...
                data_params(2,:), '\t%3.1e\n', ...
                data_params(3,:), '\t%3.1e\n', ...
                data_params(4,:), '\t%3.1e\n', ...
                data_params(5,:), '\t%3.1e\n', ...
                data_params(6,:), '\t%3.1e\n', ...
                data_params(7,:), '\t%3.1e\n', ...
                data_params(8,:), '\t%u\n', ...
                data_params(9,:), '\t%u\n', ...
                data_params(10,:), '\t%3.1e\n', ...
                data_params(11,:), '\t%u\n', ...
                data_params(12,:), '\t(%3.1e - %3.1e)\n', ...
                data_params(13,:), '\t%3.1e\n', ...
                data_params(14,:), '\t%u\n'];
    
    fprintf(fid, format_spec, wvl, wn, a0, ...
        l0, L0, delta1, deltan, N, nreals, Dz, ...
        nscr, zscr_min, zscr_max, mean(z(2:n)-z(1:n-1)),nz);
    fprintf(fid, '\nPower over slit\n');
    fprintf(fid, '%u\t', Pslit);
    
    % Close file
    fclose(fid);
end

%% SHUTDOWN COMPUTER?
wtime = 90;         % Time limit to abort shutdown

if shut
    if shut == 1
        shut_mode = 'shutdown';
        shut_cmd = 'shutdown /s';
    else
        shut_mode = 'hibernation';
        shut_cmd = 'shutdown /h';
    end
    
    fprintf('Computer %s enabled.\n', shut_mode);
    tic;
    reply_shut = 'k';
    while (~sum(strcmpi(reply_shut, {'y', 'n', ''})) && toc < wtime)
        if shut == 1
        reply_shut = waitinput('Proceed with shutdown? Y/N [Y]:    ', ...
            wtime, 's');
        else
        reply_shut = waitinput('Proceed with hibernation? Y/N [Y]:    ', ...
            wtime, 's');
        end
        fprintf('\n');
    end
    if (isempty(reply_shut) || ~sum(strcmpi(reply_shut, {'y', 'n'})))
        reply_shut = 'Y';
    end
    if strcmpi(reply_shut,'y')
        if shut == 1
            fprintf('\nShutting down.\n');
        else
            fprintf('\nHibernating.\n');
        end
        pause(3);
        shut = 0;
        system(shut_cmd);
    end
    if shut == 1
        fprintf('Shutdown aborted.\n');
    elseif shut == 2
        fprintf('Hibernation aborted.\n');
    end
end
