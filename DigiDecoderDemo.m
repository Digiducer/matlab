% DigiDecoderDemo.m
%
% This script calls the DigiDecoder function and displays the information
% returned in a dialog box
% 
clear all;
close all;

devices = DigiDecoder();

DisplayBoxStr = '';
deviceSeperatorText = sprintf('----------------------------------------------------------------------\n');
for j = 1:length(devices)
    % build up information for each device
    serialNumberText = sprintf('Serial Number: %i\n', devices(j).SN);
    calibrationDateText = sprintf('Calibration Date: %s\n', devices(j).CalDate);
    modelText = sprintf('Device: %s\n', devices(j).model);
    versionText = sprintf('Encoding version: %i\n', devices(j).version);
    if devices(j).version == 0 || devices(j).version == 1 
        sensitivityText = sprintf('Channel A Sensitivity (counts/(m/s^2)): %i\nChannel B Sensitivity (counts/(m/s^2): %i\n', devices(j).CalA, devices(j).CalB);
    elseif devices(j).version == 2 || devices(j).version == 3
        sensitivityText = sprintf('Channel A Sensitivity (counts/Volts-Pk): %i\nChannel B Sensitivity (counts/Volts-Pk): %i\n', devices(j).CalA, devices(j).CalB);
    else
        serialNumberText = sprintf('Serial Number: Not Available\n');
        sensitivityText = sprintf('Channel A Sensitivity (counts/(m/s^2)): %i\nChannel B Sensitivity (counts/(m/s^2)): %i\n', devices(j).CalA, devices(j).CalB);
    end
    message = [modelText serialNumberText versionText calibrationDateText sensitivityText deviceSeperatorText];
    DisplayBoxStr = [DisplayBoxStr message];
end

% Display Message Box with all device information
msgbox(DisplayBoxStr,'Decoded Information','modal');
