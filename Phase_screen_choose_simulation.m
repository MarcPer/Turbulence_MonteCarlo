% SUB-SCRIPT of Phase_screen_main.m
%   Request user input to choose between:
%       a) Coincidence simulation with inversion
%       b) Coincidence simulation without inversion
%       c) Intensity simulation

% Request input
choice = questdlg('Select simulation type', 'Simulation type', ...
    'Coincidences with inversion', 'Coincidences without inversion', ...
    'Intensity', 'Intensity');

% Handle response
switch choice
    case 'Coincidences with inversion'
        coinc = 1;
        invers = 1;
    case 'Coincidences without inversion'
        coinc = 1;
        invers = 0;
    case 'Intensity'
        coinc = 0;
        invers = 0;
end