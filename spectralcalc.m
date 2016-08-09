% Written by Jim Elliott for The Modal Shop, Inc.
% Modifications and documentation by Alex Lambert
function SpectrumObject = spectralcalc(timedata,offset,size,windowType)
% SpectrumObject = spectralcalc(timedata,offset,size)
% Inputs:
%   timedata: The time history of amplitudes.
%   offset: An offset, in sample counts, from the start of the timedata
%       history. Useful if you would like to skip to a later section of the
%       time data; e.g. for overlap processing. Indicates the first sample
%       used; e.g. you never have an offset of 0 since you start with the
%       1st sample at minimum.
%   size: The length of the time history, in sample counts, that you want
%       to compute spectral information over. Make sure that the size
%       accounts for your offset, if any, such that offset + size is not
%       greater than the end of your timedata.
%   window: A string specifying the window type. Valid types are:
%           'flat top','blackman-harris','hamming',and 'hann'. Feel free to
%           add your own.
%   
% Outputs:
%   SpectrumObject: A structure containing computed information. Fields
%   include:
%       Time: The time history pruned from all data as selected by inputs 
%           offset and size
%       Window: A string specifying the windowing type
%       WindowedTime: The pruned time history with windowing applied.
%       Magnitude: Ampliude corrected fft magnitudes.
%       RMS: Amplitude corrected fft RMS values.

% Compute the end sample given offset and specified length
endv = offset+size;
% Grab the time data specified from the complete history
z1=timedata(offset:endv);
% Add the selected time history to the returned SpectrumObject
SpectrumObject.Time = z1;
%NBW => noise bandwidth
if ~exist('windowType','var') ||isempty('windowType')
    windowType = 'flattop';
end
switch lower(windowType)
    case 'flattop'
        w = flattopwin(length(z1));%0.5*(1.0-cos((2.0*3.14159*(1:length(z1)))/length(z1)));
        SpectrumObject.Window = 'Flat top';
        NBW = 3.35;
    case 'flat top'
        w = flattopwin(length(z1));%0.5*(1.0-cos((2.0*3.14159*(1:length(z1)))/length(z1)));
        SpectrumObject.Window = 'Flat top';
        NBW = 3.35;
    case 'blackman-harris'
        w = blackmanharris(length(z1));
        SpectrumObject.Window = 'Blackman-Harris';
        NBW = 2.7932;
    case 'blackmanharris' 
        w = blackmanharris(length(z1));
        SpectrumObject.Window = 'Blackman-Harris';
        NBW = 2.7932;
    case 'hamming'
        w = hamming(length(z1));
        SpectrumObject.Window = 'Hamming';
        NBW = 1.36;
    case 'hamm' % does anyone call it this? just in case...
        w = hamming(length(z1));
        SpectrumObject.Window = 'Hamming';
        NBW = 1.36;
    case 'hann'
        w = hann(length(z1));
        SpectrumObject.Window = 'Hann';
        NBW=1.5;
    case 'hanning'
        w = hann(length(z1));
        SpectrumObject.Window = 'Hann';
        NBW=1.5;
    otherwise
        % default to a flat top window
        disp('Incorrect or no window specified. Defaulting to flat top.\nValid options are ''flat top'',''blackman-harris'',''hamming'',and ''hann''.');
        w = flattopwin(length(z1));%0.5*(1.0-cos((2.0*3.14159*(1:length(z1)))/length(z1)));
        SpectrumObject.Window = 'Flat top';
        NBW = 3.35;
end
%ACORR = amplitude correction
% Works for all of them
ACORR = 1/mean(w);
% Apply window
z1 = z1.*w;
SpectrumObject.WindowedTime = z1;
% Matlab FFT is 2 sided un-normalized complex fft with no correction for Window
% Amplitude. To get 1 sided mag , Take 1st half , get magnitude and
% multiply by 2/N. Then multiply by window Amplitude correction.
zf = fft(z1);
SpectrumObject.FFTdata = 2*zf(1:floor(end/2));
SpectrumObject.AmpCorr = ACORR;
SpectrumObject.NoiseBW = NBW;
% Apply magnitude corrections
SpectrumObject.Magnitude = ACORR*(abs(zf(1:floor((length(z1)/2))))*2.0)/length(z1);
SpectrumObject.RMS = sqrt((1/(2.0*NBW))*sum(SpectrumObject.Magnitude.*SpectrumObject.Magnitude));
end
