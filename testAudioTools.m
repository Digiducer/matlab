clear all

deviceReader = audioDeviceReader;
readerDevices = getAudioDevices(deviceReader)
%     hDSPAR = dsp.AudioRecorder('DeviceName', 'ASIO4ALL v2');
% hDSPAR.NumChannels = 2;
% hDSPAR.DeviceDataType = '24-bit integer';

[SN, CalA, CalB, CalDate, version, deviceNames] = DigiDecoder()
for ii=length(readerDevices):-1:1
    if length(readerDevices{ii}) < 18 || any(readerDevices{ii}(13:18) ~= '333D01')
        if length(readerDevices{ii}) >= 28 && all(readerDevices{ii}(13:28) == 'USB Audio Device')
            error('333D01 not properly detected. Unplug, close MATLAB, replug, and restart MATLAB.');
        else
            readerDevices(ii)     = [];
        end
    end
end

deviceReader = audioDeviceReader('Device',deviceNames{1},'NumChannels',2)
setup(deviceReader);
acquiredAudio = deviceReader();

plot(acquiredAudio)

