function [devices] = DigiDecoder()
%DIGIDECODER detect 'The Modal Shop' 333D01 Digiducer or
%   485B39 Digital Signal Conditioner
%   [SN, CalA, CalB, CalDate, version, ID] = DigiDecoder()
% output:
% devices - array of 333D01 and 485B39 devices found
% devices(n).SN = Device Serial Number
% devices(n).CalA = Calibration of Channel A
% devices(n).CalB = Calibration of Channel B
% devices(n).CalDate = Calibration date
% devices(n).version = device version
% devices(n).ID = Device identification

deviceReader = audioDeviceReader;
readerDevices = getAudioDevices(deviceReader);
for idx = length(readerDevices):-1:1
    if length(readerDevices{idx}) < 20 || any(readerDevices{idx}(13:18) ~= '333D01') && any(readerDevices{idx}(13:18) ~= '485B39') && any(readerDevices{idx}(13:20) ~= 'TMS 333D')
        if length(readerDevices{idx}) >= 28 && all(readerDevices{idx}(13:28) == 'USB Audio Device')
            error('333D01 or 485B39 not properly detected. Unplug, close MATLAB, replug, and restart MATLAB.');
        else
            readerDevices(idx)     = [];
        end
    end
end

% obtain sensor information
for i = 1:length(readerDevices)
    devices(i).ID = readerDevices{i};
    openParenthLocs = strfind(devices(i).ID,'(');
    closeParenthLocs = strfind(devices(i).ID,')');
    modelStr = devices(i).ID(openParenthLocs(1)+1:closeParenthLocs(1)-1);
    if modelStr(1:8) == 'TMS 333D'
        devices(i).model = 'TMS 333D';
        devices(i).CalA=33000;
        devices(i).CalB=65000;
        devices(i).version=-1;
        continue;
    end
    spaceLoc = strfind(modelStr,' ');
    devices(i).model = modelStr(1:spaceLoc-1);
    encodedStr = modelStr(spaceLoc+1:end);
    devices(i).version = str2num(encodedStr(1));
    devices(i).SN = str2num(encodedStr(2:7));
    
    switch devices(i).version
        case 0
            devices(i).CalA = str2double(encodedStr(8:12));
            devices(i).CalB = str2double(encodedStr(13:end));
        case 1
            devices(i).CalA = str2double(encodedStr(8:12));
            devices(i).CalB = str2double(encodedStr(13:17));
            dateString = encodedStr(18:end);
            yearTens = dateString(1:2);
            month = dateString(3:4);
            day = dateString(5:6);
            year = ['20',yearTens];
            devices(i).CalDate = [month,'/',day,'/',year];
        case {2, 3}
            devices(i).CalA = str2double(encodedStr(8:14));
            devices(i).CalB = str2double(encodedStr(15:21));
            if devices(i).version == 3
                devices(i).CalA = devices(i).CalA/0.05;
                devices(i).CalB = devices(i).CalB/0.05;
            end
            dateString = encodedStr(22:end);
            yearTens = dateString(1:2);
            month = dateString(3:4);
            day = dateString(5:6);
            year = ['20',yearTens];
            devices(i).CalDate = [month,'/',day,'/',year];
    end
end
if isempty(readerDevices)
    error('No 333D01 or 485B39 was detected. Unplug, close MATLAB, replug, and restart MATLAB.');
end