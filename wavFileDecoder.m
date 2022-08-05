function [ CalDate, SN, CalA, CalB, returnVal ] = wavFileDecoder( Wav_file )
%wavFileDecoder - This file reads the format of a wav_file and returns back
% the calibration information
%
%   The program first confirms it's a wave file, finds the data chunk for
%   calibration, extracts that data, and sends it to be plotted and
%   confirmed

%open the file to read byte by byte
Sensitivity_file = fopen(Wav_file, 'r');

%iterate through and grab the header data
for header_count = 1:4
    
    array = fread(Sensitivity_file, 4, '*char');
    if header_count == 1
        
        RIFF_comp = strcat(array(1), array(2), array(3), array(4));
    elseif header_count == 2
        
        Arbitrary_size = strcat(array(1), array(2), array(3), array(4));
    elseif header_count == 3
        
        WAVE_comp = strcat(array(1), array(2), array(3), array(4));
    else
        
        FMT_comp = strcat(array(1), array(2), array(3), array(4));
    end
    
    header_count = header_count + 1;
end
%move on only if the correct wav file header was found
if strcmp(RIFF_comp, 'RIFF') == 1 && strcmp(WAVE_comp, 'WAVE') == 1 && strcmp(FMT_comp, 'fmt') == 1
    bytes = fread(Sensitivity_file, 1, '*int');
    
    if bytes == 0
        errordlg('There is no remaining data to read or seek for. ABORT');
        returnVal = -1;
        return
    end
    
    status = fseek(Sensitivity_file, bytes, 'cof');
   
    if status == -1
        errordlg('There is no remaining data to read or seek for. ABORT');
        returnVal = -1;
        return
    end
    
    array = fread(Sensitivity_file, 4, '*char');
    
    array = strcat(array(1), array(2), array(3), array(4));
    %continue to look for Cal1 data
    while  strcmp(array, 'CAL1' ) ~= 1
        
        bytes = fread(Sensitivity_file, 1, '*int');
        
        fseek(Sensitivity_file, bytes, 'cof');
        
        array = fread(Sensitivity_file, 4, '*char');
        %if there is no calibration data, set these as defaults and return
        %0
        if length(array) ~= 4
           returnVal = 0;
           SN = 'Not Available';
           %input nominal values
           CalA = 33000;
           CalB = 65000;
           CalDate = 'Not Available';
           return
        end
        
        array = strcat(array(1), array(2), array(3), array(4));

    end
    %%%%% GUYS GUYS! WE FOUND CAL1!!!
    bytes = fread(Sensitivity_file, 1, '*int');
    
    Cal1_content = transpose(fread(Sensitivity_file, bytes, '*char'));
    
    Cal1_content = strsplit(Cal1_content, ' ');
    
    if strcmp(Cal1_content(1), '333D01') ~= 1
        
        errordlg('This is not the 333D01 sensor. ABORT!');
        
    end
    
    encodedStr = cell2mat(Cal1_content(3));
    %make sure it's the correct version
    if str2num(encodedStr(1)) ~= 1
        
        errordlg('This is not formatted to version 1. ABORT!');
    end
    %place all the data in arrays and return a 1
    SN = str2num(encodedStr(2:7));
    
    CalA = str2num(encodedStr(8:12));
    CalB = str2num(encodedStr(13:17));
    dateString = encodedStr(18:end);
    yearTens = dateString(1:2);
    month = dateString(3:4);
    day = dateString(5:6);
    year = ['20',yearTens];
    CalDate = [month,'/',day,'/',year];

    returnVal = 1;
else
    %if not formatted correctly, return -1
    errordlg('This is not a properly formatted wav file. ABORT!');
    returnVal = -1;
end


