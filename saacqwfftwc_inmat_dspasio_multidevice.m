% Given a sample rate, block size, and average count, captures data using
% the DSP toolkit hook to ASIO and computes a fourier spectrum.
% Capable of capturing from multiple devices at the same time. However, in
% order to prevent a buffer overrun, capture is paused while spectrum
% calculations are made. This means that the displayed infomation is
% real-time and multiple device/channel, but not continuous.
% Continuous capture is possible, for 1-2 devices.
clear y yHist xt1 xt2 xs1 xs2 xavg1 xavg2 xtaxis xfaxis
if exist('deviceReader','var')
%     release(hDSPAR);
%     clear hDSPAR
    release(deviceReader);
end
prompt = {'Averages',sprintf('Sample rate - Valid options are\n8000,11025,16000,22050,32000,44100,48000'),'Block size (Max 65536)','Number of Devices','Sensitivity (Ch. A),mV/g or AD counts / (m/s^2)','Sensitivity (Ch. B), mV/g or AD counts / (m/s^2)','Window'};
% Important note: The 333D01 only supports sample rates of:
% 8000Hz, 11025Hz, 16000Hz, 22050Hz, 32000Hz, 44100Hz, 48000Hz
answer  = inputdlg(prompt,'Acquisition',1,{'30','48000','16384','1','36738','71964','flattop'});
ansarr = cellfun(@str2num,answer(1:end-1));
window = answer{end};
averages = ansarr(1);
sampleRates = ansarr(2);
blockSizes = ansarr(3);
numDevs = ansarr(4);
sensitivityA_user = ansarr(5);
sensitivityB_user = ansarr(6);
% set this to true to accumulate all time data
plotAllTime = false;
% set this to true to never stop filling the buffer. 
% False means data collection stops while spectrums are calculated.
continuousMode = true;
if sensitivityA_user > 10000
    %Entered counts / m/s^2
    sensitivity_ctpmps2A = sensitivityA_user;
elseif sensitivityA_user > 40 && sensitivityA_user < 80
    % is incorrect old unit, mV/m/s^2
    sensitivity_ctpmps2A = sensitivityA_user*2^23/(1000*sqrt(2)*9.80665);
end
if sensitivityB_user > 10000
    %Entered counts / m/s^2
    sensitivity_ctpmps2B = sensitivityB_user;
elseif sensitivityB_user > 100 && sensitivityB_user < 140
    % is incorrect old unit, mV/m/s^2
    sensitivity_ctpmps2B = sensitivityB_user*2^23/(1000*sqrt(2)*9.80665);
end
scaleFactorA = 2^23 / sensitivity_ctpmps2A;
scaleFactorB = 2^23 / sensitivity_ctpmps2B;
% initialize recording device
% hDSPAR = dsp.AudioRecorder('DeviceName','ASIO4ALL v2');
% hDSPAR.NumChannels = 2;
% hDSPAR.DeviceDataType ='24-bit integer';
% If buffer overrun occurs, you can try increasing the buffer size or the
% queue duration. However, the essential problem is that MATLAB takes too
% long between buffer dumps to variable (caling step(h), now deviceReader()), and too much data
% accumulates. Moving your block size closer to or higher than the sample
% rate can help as well; e.m/s^2. a sample rate of 48kHz and a block of 8192
% samples is more likely to cause overrun than 8000 Hz and 8192 samples.
% You can increase the queue duration until your memory fills up if you
% need the higher sample rate. Doing so will greatly increase latency,
% however. Check documentation for dsp.AudioRecorder for more information
% on latency.
% hDSPAR.BufferSizeSource ='Auto';
% h.BufferSizeSource ='Property';
% h.BufferSize = 4096;
% hDSPAR.QueueDuration = 15;
% hDSPAR.SampleRate = sampleRates;
% hDSPAR.SamplesPerFrame = blockSizes;
deviceReader = audioDeviceReader(sampleRates,blockSizes,'Device','Default',...
                                 'BitDepth','24-bit integer','NumChannels',2);
% We need to give ASIO4ALL a chance to get set up for our planned
% usage.
% To do so, click on the green arrow icon in the taskbar, and open the
% advanced menu. Click in the blue square icons to enable or disable a
% device. Disable all windows devices, and enable all of the 333D01
% USB Smart Sensors that are connected that you would like to use.
% If a devices is shown in an unusual state, try powering the device
% on/off, reconnecting the device, and then restarting MATLAB.

% Give ASIO4ALL a second to get initialized
% grab some data to start up the icon which gives user access to
% settings
% release(hDSPAR);
% [y] = step(hDSPAR);
release(deviceReader);
y=deviceReader();
box = msgbox('Please select the ASIO4ALL icon and configure your devices for use.','Configure ASIO4ALL');
% wait for user action to continue
uiwait(box);
% release control of the device until we begin data capture
% release(hDSPAR);
release(deviceReader);
pause(1);
hDSPAR.NumChannels = 2*numDevs;
% initialize these
avgFreqs = zeros(length(sampleRates),1);
peakAvgMags = zeros(length(sampleRates),1);
for devIdx = 1:numDevs
    figh(devIdx) = figure;
    % set up handles to subplot regions
    sp(devIdx,1) = subplot(3,1,1);     % Time
    sp(devIdx,2) = subplot(3,1,2);     % Frequency
    sp(devIdx,3) = subplot(3,1,3);     % Frequency avg
    % Axis scaling
    % time axis - divisions by 1/sampleRate
%     xtaxis(devIdx,:) = (1/hDSPAR.SampleRate)*(0:hDSPAR.SamplesPerFrame-1);
    xtaxis(devIdx,:) = (1/sampleRates)*(0:blockSizes-1);
    % frequency axis - divisions follow frequency resolution defined by
    % sampleRate / blockSize
%     xfaxis(devIdx,:) = (hDSPAR.SampleRate/(hDSPAR.SamplesPerFrame))*(0:(hDSPAR.SamplesPerFrame/2)-1);
    xfaxis(devIdx,:) = (sampleRates/(blockSizes))*(0:(blockSizes/2)-1);
end
if plotAllTime
    % initialize the all-time-history variable for speed
    yHist = zeros(blockSizes*averages,2*numDevs);
end
% release(hDSPAR); % make sure devices are not already in use before starting capture loop
release(deviceReader);
for i=1:averages
    % grab data from device, start buffer loading
%     y = step(hDSPAR);
    y=deviceReader();
    if plotAllTime
        % copy data to all-time buffer.
        yHist((i-1)*blockSizes+1:(i*blockSizes),:) = y(:,:);
    end
    if ~continuousMode
%         release(hDSPAR); % stop data collection into the buffer to allow time for
                    % spectral calculations to complete
        release(deviceReader);
    end
    for devIdx = 1:numDevs
        % switch to a new figure window for each device
        set(0,'CurrentFigure',figh(devIdx));
        % Select time data display
        subplot(sp(devIdx,1));
        % Time history plot
        plot(xtaxis,y(:,1+(devIdx-1)*2)*scaleFactorA,xtaxis,y(:,2+(devIdx-1)*2)*scaleFactorB);
        xlabel('Time (s)');
        ylabel(sprintf('Acceleration\n(m/s^2)'));
		legend('Ch 1','Ch 2');
        xlim([0 max(xtaxis(devIdx,:))]);
        % Compute amplitude of time history for this average sample
        timeDomMax1(devIdx,i) = max(y(:,1+(devIdx-1)*2)*scaleFactorA);
        timeDomMax2(devIdx,i) = max(y(:,2+(devIdx-1)*2)*scaleFactorB);
        %Select frequency display
        subplot(sp(devIdx,2));
        % get time history for channel A
        xt1(devIdx,:) = y(:,1+(devIdx-1)*2);
        % Compute spectrum for channel A
%         xs1(devIdx,:)=spectralcalc(xt1(devIdx,:)',1,hDSPAR.SamplesPerFrame-1,window); % scaling the halved channel to get the correct vibration amplitude
        xs1(devIdx,:)=spectralcalc(xt1(devIdx,:)',1,blockSizes-1,window); % scaling the halved channel to get the correct vibration amplitude
        % get time history for channel B
        xt2(devIdx,:) = y(:,2+(devIdx-1)*2);
        % Compute spectrum for channel B
%         xs1(devIdx,:)=spectralcalc(xt1(devIdx,:)',1,hDSPAR.SamplesPerFrame-1,window); % scaling the halved channel to get the correct vibration amplitude
        xs2(devIdx,:)=spectralcalc(xt2(devIdx,:)',1,blockSizes-1,window); % scaling the halved channel to get the correct vibration amplitude
        % averaging
        if i == 1
            % Channel A average for first sample is itself
            xavg1(devIdx,:) = xs1(devIdx,:).Magnitude;
            % Channel B average for first sample is itself
            xavg2(devIdx,:) = xs2(devIdx,:).Magnitude;
        else
            % Channel A average
            xavg1(devIdx,:) = xavg1(devIdx,:)+(xs1(devIdx,:).Magnitude)';
            % Channel B average
            xavg2(devIdx,:) = xavg2(devIdx,:)+(xs2(devIdx,:).Magnitude)';
        end
        % Compute peak statistics to adjust display
        [~,peakFreqInd(devIdx)] = max(xs2(devIdx,1).Magnitude*scaleFactorB);
        peakFreq(devIdx) = xfaxis(devIdx,peakFreqInd(devIdx)); 
        peakFreqMag(devIdx) = xs2(devIdx,1).Magnitude(peakFreqInd(devIdx))*scaleFactorB;
        semilogx(xfaxis(devIdx,:),20.*log10(xs1(devIdx,1).Magnitude*scaleFactorA),xfaxis(devIdx,:),20.*log10(xs2(devIdx,1).Magnitude*scaleFactorB));
        xlabel('Frequency (Hz)'),ylabel(sprintf('Acceleration\n(dB of m/s^2)'));
		% Keep legend out of the way
        if peakFreq < 500
            legend('Ch 1','Ch 2','location','NorthEast');
        elseif peakFreq >= 500
            legend('Ch 1','Ch 2','location','NorthWest');
        end;
        grid on;
        xlim([0 10000]);
        ylim([-70 20*log10(peakFreqMag(devIdx))+10]); 
		
        % averaged spectrum       
        subplot(sp(devIdx,3));
        semilogx(xfaxis(devIdx,:),20.*log10(xavg1(devIdx,:)*scaleFactorA/i),xfaxis(devIdx,:),20.*log10(xavg2(devIdx,:)*scaleFactorB/i));
        xlabel('Frequency (Hz)'),ylabel(sprintf('Avg Acceleration\n(dB of m/s^2)\nAvg #: %i',i));
        % keep legend out of the way
        if peakFreq < 500
            legend('Ch 1','Ch 2','location','NorthEast');
        elseif peakFreq >= 500
            legend('Ch 1','Ch 2','location','NorthWest');
        end;
        grid on;
        xlim([0 10000]);
        ylim([-70 20*log10(peakFreqMag(devIdx))+10]);
        % force draw
        drawnow();
    end
end
for devIdx = 1:numDevs
    % Compute averaged statistics on identified peak frequencies and their
    % magnitudes, as well as time history amplitude
    [~,avgPeakFreqInd1(devIdx)] = max(xavg1(devIdx,:)*scaleFactorA/averages);
    [~,avgPeakFreqInd2(devIdx)] = max(xavg2(devIdx,:)*scaleFactorB/averages);
    peakAvgFreq1(devIdx) = xfaxis(devIdx,avgPeakFreqInd1(devIdx));
    peakAvgFreq2(devIdx) = xfaxis(devIdx,avgPeakFreqInd2(devIdx));
    avgPeakAmplitude1(devIdx) = 20*log10(max(xavg1(devIdx,:)*scaleFactorA)/averages);
    avgPeakAmplitude2(devIdx) = 20*log10(max(xavg2(devIdx,:)*scaleFactorB)/averages);
    avgTimePeakAmplitude1(devIdx) = mean(timeDomMax1(devIdx));
    avgTimePeakAmplitude2(devIdx) = mean(timeDomMax2(devIdx));
    % print out summary
    fprintf('Ch.A Peak Frequency: %.2f Hz. Ch.A Peak magnitude: %.3f dB acceleration. Ch.A Peak time amp: %.4f m/s^2\n',peakAvgFreq1(devIdx),avgPeakAmplitude1(devIdx),avgTimePeakAmplitude1(devIdx));
    fprintf('Ch.B Peak Frequency: %.2f Hz. Ch.B Peak magnitude: %.3f dB acceleration. Ch.B Peak time amp: %.4f m/s^2\n\n',peakAvgFreq2(devIdx),avgPeakAmplitude2(devIdx),avgTimePeakAmplitude2(devIdx));
end
% create a plot of the total time history signal for all devices, separated
% into channel As and channel Bs
if plotAllTime
    figure(numDevs+1);
    cla;
    figure(numDevs+2);
    cla;
    for idx = 1:numDevs
        %create new figure for all devices, 
        % one for ch a and another for ch b
        figure(numDevs +1);
        hold all;
%         plot(1/hDSPAR.SampleRate.*[1:length(yHist(:,1))],yHist(:,2*(idx)-1)*scaleFactorA);
        plot(1/sampleRates.*[1:length(yHist(:,1))],yHist(:,2*(idx)-1)*scaleFactorA);
        title('Channel A');
        xlabel('Time (s)');
        ylabel('Acceleration (m/s^2)');
        figure(numDevs +2);
        hold all;
%         plot(1/hDSPAR.SampleRate.*[1:length(yHist(:,1))],yHist(:,2*idx)*scaleFactorB);
        plot(1/sampleRates.*[1:length(yHist(:,1))],yHist(:,2*idx)*scaleFactorB);
        title('Channel B');
        xlabel('Time (s)');
        ylabel('Acceleration (m/s^2)');
    end
end
% remove unimportant variables to avoid clogging workspace
clear prompt answer ansarr plotAllTime