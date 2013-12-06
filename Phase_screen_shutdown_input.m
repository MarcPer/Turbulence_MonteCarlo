% SUB-SCRIPT of Phase_screen_main.m
%   Request user input when computer is scheduled for shutdown or
%   hibernate.
%   Choices are between:
%       0: DEACTIVATE SHUTDOWN (OR HIBERNATE)
%       1: PROCEED WITH SHUTDOWN (OR HIBERNATE)

% Input error check
arg_defined = exist('shut', 'var'); % Was SHUT defined?
if ~arg_defined
    error('PhaseScreen_shut:noarg', 'Variable SHUT must be defined.');
end
if ~isscalar(shut)
    error('PhaseScreen_shut:arg', ['Variable SHUT must be either', ...
        ' 0, 1 or 2.']);
end


if shut         % Shutdown was activated - confirm with user
    if all( ~(shut == [0; 1; 2]))   % Check if argument is correct
        error('PhaseScreen_shut:arg', ['Variable SHUT must be either', ...
            ' 0, 1 or 2.']);
    end
    
    % Shutdown or hibernate?
    if shut == 1
        pmpt_shut = ['Computer is scheduled for shutdown. Press', ...
            ' CANCEL to deactivate shutdown.'];
    end
    if shut == 2
        pmpt_shut = ['Computer is scheduled for hibernation. Press', ...
            ' CANCEL to deactivate hibernation.'];
    end
    
    % Request input
    choice_shut = questdlg(pmpt_shut, 'Confirm Shutdown/Hibernation', ...
        'OK', 'Cancel', 'Cancel');

    % Handle response
    switch choice_shut
        case 'Cancel'
            shut = 0;
        case 'OK'
            if shut == 1
                fprintf('Shutdown confirmed.\n');
            else
                fprintf('Hibernation confirmed.\n');
            end
    end
end