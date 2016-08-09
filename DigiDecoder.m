function [SN, CalA, CalB, CalDate, version] = DigiDecoder()
% close all;
a = audiodevinfo;
allDevs = {a.input.Name};
for idx = length(allDevs):-1:1
    if allDevs{idx}(13:18) ~= '333D01'
        if allDevs{idx}(13:28) == 'USB Audio Device'
            error('333D01 not properly detected. Unplug, close MATLAB, replug, and restart MATLAB.');
        else
            allDevs(idx)     = [];
        end
    end    
end
%[selections, ok] = listdlg('PromptString',sprintf('Select which device(s) to decode. (Ctrl+click or click and drag for multiple select)'),'SelectionMode','multiple','ListString',allDevs,'ListSize',[400 100]);
SN = 0;
CalA = 0;
CalB = 0;
CalDate = '';
    for i = 1:length(allDevs)
    inputStr = allDevs{i};
    totalStrLen = length(inputStr);
    openParenthLocs = findstr(inputStr,'(');
    closeParenthLocs = findstr(inputStr,')');
    interface = inputStr(openParenthLocs(2)+1:closeParenthLocs(2)-1);
    modelStr = inputStr(openParenthLocs(1)+1:closeParenthLocs(1)-1);
    spaceLoc = findstr(modelStr,' ');
    model = modelStr(1:spaceLoc-1);
    encodedStr = modelStr(spaceLoc+1:end);
    version(i) = str2num(encodedStr(1));
    
    switch version(i)
        case 0
            SN(i) = str2num(encodedStr(2:7));
            CalA(i) = str2num(encodedStr(8:12));
            CalB(i) = str2num(encodedStr(13:end));
        case 1
            SN(i) = str2num(encodedStr(2:7));
            CalA(i) = str2num(encodedStr(8:12));
            CalB(i) = str2num(encodedStr(13:17));
            dateString = encodedStr(18:end);
            yearTens = dateString(1:2);
            month = dateString(3:4);
            day = dateString(5:6);
            year = ['20',yearTens];
            CalDate{i} = [month,'/',day,'/',year];
    end
    end
if length(allDevs) > 1
    msgbox(sprintf('Warning, multiple 333D01s detected. Defaulting to serial number %i',SN(1)),'More than one device');
end