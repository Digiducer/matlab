% Given a sample rate, block size, and average count, captures data using
% the DSP toolkit hook to ASIO and computes a fourier spectrum.
clear y yHist xt1 xt2 xs1 xs2 xavg1 xavg2 xtaxis xfaxis
if exist('deviceReader','var')
    % make sure to clear h after switching from other capture methods
%     release(hDSPAR);
%     clear hDSPAR
    release(deviceReader);
end
prompt = {'Averages',sprintf('Sample rate - Valid options are\n8000,11025,16000,22050,32000,44100,48000'),'Block size (Max 65536)','Window'};
% Important note: The 333D01 only supports sample rates of:
% 8000Hz, 11025Hz, 16000Hz, 22050Hz, 32000Hz, 44100Hz, 48000Hz
answer  = inputdlg(prompt,'Acquisition',1,{'30','48000','16384','flattop'});
ansarr = cellfun(@str2num,answer(1:end-1));
window = answer{end};
NUM_AVERAGES = ansarr(1);
sampleRates = ansarr(2);
blockSizes = ansarr(3);
[SN, calA, calB, calDate, vers] = DigiDecoder;
sensitivityA_user = calA(1);
sensitivityB_user = calB(1);
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
% set this to true to accumulate all time data
plotAllTime = false;
% set this to true to never stop filling the buffer. 
% False means data collection stops while spectrums are calculated.
continuousMode = true;
% initialize recording device
% hDSPAR = dsp.AudioRecorder('DeviceName','Default');
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
% hDSPAR.BufferSizeSource ='Property';
% hDSPAR.BufferSize = 8192;
% hDSPAR.QueueDuration = 5;
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
% [y] = step(hDSPAR);
y=deviceReader();
box = msgbox('Please select the ASIO4ALL icon and configure your devices for use.','Configure ASIO4ALL');
% wait for user action to continue
uiwait(box);
% release control of the device until we begin data capture
% release(hDSPAR);
release(deviceReader);
pause(1);
% initialize these
avgFreqs = zeros(length(sampleRates),1);
peakAvgMags = zeros(length(sampleRates),1);

figh = figure;
% set up handles to subplot regions
sp(1) = subplot(3,1,1);     % Time
sp(2) = subplot(3,1,2);     % Frequency
sp(3) = subplot(3,1,3);     % Frequency avg
% Axis scaling
% time axis - divisions by 1/sampleRate
% xtaxis = (1/hDSPAR.SampleRate)*(0:hDSPAR.SamplesPerFrame-1);
xtaxis = (1/sampleRates)*(0:blockSizes-1);
% frequency axis - divisions follow frequency resolution defined by
    % sampleRate / blockSize
% xfaxis = (hDSPAR.SampleRate/(hDSPAR.SamplesPerFrame))*(0:(hDSPAR.SamplesPerFrame/2)-1);
xfaxis = (sampleRates/(blockSizes))*(0:(blockSizes/2)-1);
for i=1:NUM_AVERAGES
	% grab data from device, start buffer loading
% 	y = step(hDSPAR);
    y=deviceReader();
	if plotAllTime
        % copy data to all-time buffer.
        yHist((i-1)*blockSizes+1:(i*blockSizes),:) = y(:,:);
    end
    if ~continuousMode
%         release(hDSPAR); % stop data collection into the buffer to allow time for
        release(deviceReader);
                    % spectral calculations to complete
    end
	set(0,'CurrentFigure',figh);
	% Select time data display
    if i == 1
        subplot(sp(1));
        p1 = plot(xtaxis,y(:,1).*scaleFactorA,xtaxis,y(:,2).*scaleFactorB);
        xlabel('Time (s)');
        ylabel(sprintf('Acceleration\n(m/s^2)'));
        xlim([0 max(xtaxis)]);
    else
        set(p1(1),'XData',xtaxis,'YData',y(:,1).*scaleFactorA);
        set(p1(2),'XData',xtaxis,'YData',y(:,2).*scaleFactorB);
    end
	% Compute amplitude of time history for this average sample
	timeDomMax1(i) = max(y(:,1)*scaleFactorA);
	timeDomMax2(i) = max(y(:,2)*scaleFactorB);
	% get time history for channel A
	xt1 = y(:,1);
	% Compute spectrum for channel A
% 	xs1=spectralcalc(xt1,1,hDSPAR.SamplesPerFrame-1,window); % scaling the halved channel to get the correct vibration amplitude
	xs1=spectralcalc(xt1,1,blockSizes-1,window); % scaling the halved channel to get the correct vibration amplitude
	% get time history for channel B
	xt2 = y(:,2);
	% Compute spectrum for channel B
% 	xs2=spectralcalc(xt2,1,hDSPAR.SamplesPerFrame-1,window);
	xs2=spectralcalc(xt2,1,blockSizes-1,window);
	% averaging
	if i == 1
		% Channel A average for first sample is itself
		xavgsum1 = xs1.Magnitude;
		% Channel B average for first sample is itself
		xavgsum2 = xs2.Magnitude;
	else
		% Channel A average
		xavgsum1 = xavgsum1+xs1.Magnitude;
		% Channel B average
		xavgsum2 = xavgsum2+xs2.Magnitude;
	end
	% Compute peak statistics to adjust display
   [~,peakFreqInd] = max(xs2(1).Magnitude*scaleFactorB);
    peakFreq = xfaxis(peakFreqInd); 
    peakFreqMag = xs2(1).Magnitude(peakFreqInd)*scaleFactorB;
   if i == 1
        subplot(sp(2));
        p2 = semilogx(xfaxis,20.*log10(xs1(1).Magnitude*scaleFactorA),xfaxis,20.*log10(xs2(1).Magnitude*scaleFactorB));
        xlabel('Frequency (Hz)'),ylabel(sprintf('Acceleration\n(dB of m/s^2)'));
        % Keep legend out of the way
        if peakFreq < 500
            legend('Ch 1','Ch 2','location','NorthEast');
        elseif peakFreq >= 500
            legend('Ch 1','Ch 2','location','NorthWest');
        end;
        grid on;
        xlim([0 10000]);
        ylim([-70 mag2db(peakFreqMag)+10]); 
    else
        set(p2(1),'XData',xfaxis,'YData',20.*log10(xs1(1).Magnitude*scaleFactorA));
        set(p2(2),'XData',xfaxis,'YData',20.*log10(xs2(1).Magnitude*scaleFactorB));
    end
	
    if i == 1
        subplot(sp(3));
        p3 = semilogx(xfaxis,20.*log10(xavgsum1/min(i,NUM_AVERAGES)*scaleFactorA),xfaxis,20.*log10(xavgsum2/min(i,NUM_AVERAGES)*scaleFactorB));
        xlabel('Frequency (Hz)'),ylabel(sprintf('Avg Acceleration\n(dB of m/s^2)\nAvg #: %i',i));
        % keep legend out of the way
        if peakFreq < 500
            legend('Ch 1','Ch 2','location','NorthEast');
        elseif peakFreq >= 500
            legend('Ch 1','Ch 2','location','NorthWest');
        end;
        grid on;
        xlim([0 10000]);
        ylim([-70 mag2db(peakFreqMag)+10]); 
    else
        set(p3(1),'XData',xfaxis,'YData',20.*log10(xavgsum1/min(i,NUM_AVERAGES)*scaleFactorA));
        set(p3(2),'XData',xfaxis,'YData',20.*log10(xavgsum2/min(i,NUM_AVERAGES)*scaleFactorB));
        ylabel(sprintf('Avg Acceleration\n(dB of m/s^2)\nAvg #: %i',min(i,NUM_AVERAGES)));
    end
    % force draw
    drawnow();
end
% Compute averaged statistics on identified peak frequencies and their
    % magnitudes, as well as time history amplitude
    [~,avgPeakFreqInd1] = max(xavgsum1*scaleFactorA/NUM_AVERAGES);
    [~,avgPeakFreqInd2] = max(xavgsum2*scaleFactorB/NUM_AVERAGES);
    peakAvgFreq1 = xfaxis(avgPeakFreqInd1);
    peakAvgFreq2 = xfaxis(avgPeakFreqInd2);
    avgPeakAmplitude1 = 20*log10(max(xavgsum1)*scaleFactorA/NUM_AVERAGES);
    avgPeakAmplitude2 = 20*log10(max(xavgsum2)*scaleFactorB/NUM_AVERAGES);
    avgTimePeakAmplitude1 = mean(timeDomMax1);
    avgTimePeakAmplitude2 = mean(timeDomMax2);
    % print out summary
    fprintf('Ch.A Peak Frequency: %.2f Hz. Ch.A Peak magnitude: %.3f dB acceleration. Ch.A Peak time amp: %.4f m/s^2\n',peakAvgFreq1,avgPeakAmplitude1,avgTimePeakAmplitude1);
    fprintf('Ch.B Peak Frequency: %.2f Hz. Ch.B Peak magnitude: %.3f dB acceleration. Ch.B Peak time amp: %.4f m/s^2\n\n',peakAvgFreq2,avgPeakAmplitude2,avgTimePeakAmplitude2);

	% create a plot of the total time history signal for all devices, separated
% into channel As and channel Bs
if plotAllTime
    figure(2);
    cla;
    figure(3);
    cla;
    %create new figure for all devices, 
    % one for ch a and another for ch b
    figure(2);
    hold all;
    plot(1/hDSPAR.SampleRate.*[1:length(yHist(:,1))],yHist(:,1)*scaleFactorA);
    title('Channel A');
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
    figure(3);
    hold all;
    plot(1/hDSPAR.SampleRate.*[1:length(yHist(:,2))],yHist(:,2)*scaleFactorB);
    title('Channel B');
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
end
% remove unimportant variables to avoid clogging workspace
clear prompt answer ansarr plotAllTime