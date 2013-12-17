classdef UserInput<handle
    properties
        enumShutdown;
        isFourthOrder;
        isInverted;
    end
    
    methods(Access = public)
       function getSimulationType(obj)
            choice = questdlg('Select simulation type', 'Simulation type', ...
                'Coincidences with inversion', 'Coincidences without inversion', ...
                'Intensity', 'Intensity');
            
            switch choice
                case 'Coincidences with inversion'
                    obj.isFourthOrder = 1;
                    obj.isInverted = 1;
                case 'Coincidences without inversion'
                    obj.isFourthOrder = 1;
                    obj.isInverted = 0;
                case 'Intensity'
                    obj.isFourthOrder = 0;
                    obj.isInverted = 0;
            end
        end 
        
    end
    
    methods(Static)
        function isAbort = abortWhenConstraintFail(fail, shutCtrl)
            isAbort = 0;
            if ~fail
                fprintf('Press any key to continue...\n');
                pause;
                return
            end
            switch( questdlg('Abort?','Constraint analysis failed', 'Yes', 'No', 'Yes'));
                case 'Yes'
                    isAbort = 1;
            end
            UserInput.shutDownAborted(isAbort, shutCtrl);
        end
        function shutDownAborted(isAbort, shutCtrl)
            if isAbort
                switch shutCtrl
                    case 1
                        fprintf('Shutdown aborted...\n');
                    case 2
                        fprintf('Hibernation aborted...\n');
                end
            end
        end
        function h = createWaitBar()
            h = waitbar(0, '0%', 'Name', 'Simulating turbulence...', ...
                'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0);
        end
        function updateWaitBar(waitBarHandle,iRe, nRe)
            waitbar(iRe/nRe, waitBarHandle, sprintf('%3.1f%%',iRe/nRe*100)); 
        end
        function abort = isAborted(waitBarHandle)
            abort = 0;
            if getappdata(waitBarHandle,'canceling')
                abort = 1;
            end
        end
        function printOutProgress(statusString, idx, fIdx)
            fprintf('%s: %u of %u\n', statusString, idx, fIdx);
        end
    end
    
    methods
        function set.enumShutdown(obj,value)
            if all( ~(value == [0; 1; 2]))   % Check if argument is correct
                error('PhaseScreen_shut:arg', ['Variable SHUT must be either', ...
                    ' 0, 1 or 2.']);
            end
            
            obj.enumShutdown = 0;
            if (value == 1)
                shutMode = 'SHUTDOWN';
            end
            if (value == 2)
                shutMode = 'HIBERNATION';
            end
            
            if value
                prompt = sprintf('Computer is scheduled for %s. Press CANCEL to unschedule.\n', shutMode);
                choiceShut = questdlg(prompt, 'Confirm Shutdown/Hibernation', ...
                    'OK', 'Cancel', 'Cancel');
            end
            
            switch choiceShut
                case 'Cancel'
                    obj.enumShutdown = 0;
                    fprintf('%s cancelled.\n', shutMode);
                case 'OK'
                    fprintf('%s confirmed.\n', shutMode);
            end
            
        end
    end
    
end