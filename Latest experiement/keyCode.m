function keyNumCode = keyCode( keyName, context )
% function keyNumCode = keyCode( keyName, context )
% Lookup of the numerical code for a KEY_NAME string in a particular
% context, e.g. the NSPN Dell Latitudes or a scanner ...

%% Default arguments:
if nargin < 2
    context = 1; 
end;

%% contexts to be catered for :
Dell_Latitude = 1;
% Add lines about a particular device (keypad) here ...

%% Lookup:
if  context == Dell_Latitude
    switch keyName
        case 'RETURN'
            keyNumCode = 59;
        case 'SPACE'
            keyNumCode = 71;
        case 'UP_ARROW'   
            keyNumCode = 95;
        case 'LEFT_ARROW'   
            keyNumCode = 97;
        case 'RIGHT_ARROW'
            keyNumCode = 98;
        case 'DOWN_ARROW'
            keyNumCode = 100;
        case 'ESCAPE'
            keyNumCode = 52;
        case 'F'    % Used sometimes for 'move forward' or 'finish this part'
            keyNumCode = 6;
        case 'G'    % Used in MVS3
            keyNumCode = 7;
        case 'S'
            keyNumCode = 19;
        case '0'
            keyNumCode = 27;
        case '1'
            keyNumCode = 28;
        case '2'
            keyNumCode = 29;
        case '3'
            keyNumCode = 30;
        case '4'
            keyNumCode = 31;
        case '5'
            keyNumCode = 32;
        case '6'
            keyNumCode = 33;
        case '7'
            keyNumCode = 34;
        case '8'
            keyNumCode = 35;
        case '9'
            keyNumCode = 36;
        otherwise
            error('Valid names are: ''ESCAPE'',''RETURN'',''SPACE'',''UP_ARROW'',''LEFT_ARROW'',''RIGHT_ARROW'',''DOWN_ARROW'',''F'',''G'',''S'',''0-9''');
    end % of choices available for context Dell_Latitude
else
    error('function keyCode(keyName, context) unprepared for the context provided');
end
    
return; % of whole function