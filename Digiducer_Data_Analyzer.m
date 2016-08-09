function varargout = Digiducer_Data_Analyzer(varargin)
% DIGIDUCER_DATA_ANALYZER MATLAB code for Digiducer_Data_Analyzer.fig
%      DIGIDUCER_DATA_ANALYZER, by itself, creates a new DIGIDUCER_DATA_ANALYZER or raises the existing
%      singleton*.
%
%      H = DIGIDUCER_DATA_ANALYZER returns the handle to a new DIGIDUCER_DATA_ANALYZER or the handle to
%      the existing singleton*.
%
%      DIGIDUCER_DATA_ANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGIDUCER_DATA_ANALYZER.M with the given input arguments.
%
%      DIGIDUCER_DATA_ANALYZER('Property','Value',...) creates a new DIGIDUCER_DATA_ANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Digiducer_Data_Analyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Digiducer_Data_Analyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Digiducer_Data_Analyzer

% Last Modified by GUIDE v2.5 23-Feb-2015 10:59:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Digiducer_Data_Analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @Digiducer_Data_Analyzer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Digiducer_Data_Analyzer is made visible.
function Digiducer_Data_Analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Digiducer_Data_Analyzer (see VARARGIN)

% Choose default command line output for Digiducer_Data_Analyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Digiducer_Data_Analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Digiducer_Data_Analyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% % --- Executes on selection change in blockSize.
function blockSize_Callback(hObject, eventdata, handles)
% % hObject    handle to blockSize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns blockSize contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from blockSize
% blockSize = 2^(5+get(handles.blockSize, 'Value'));

% --- Executes during object creation, after setting all properties.
function blockSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set up options for block size
set(hObject, 'String', {'64', '128', '256', '512', '1024', '2048', '4096', '8192', '16384', '32768', '65536'});
%initialize block size to 1024
set(hObject, 'Value', 5);

% --- Executes on selection change in windowType.
function windowType_Callback(hObject, eventdata, handles)
% hObject    handle to windowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns windowType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from windowType

% --- Executes during object creation, after setting all properties.
function windowType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set up options for window type
set(hObject, 'String', {'flattop', 'hanning', 'hann', 'hamming', 'blackman-harris'});

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%resets and clears all output boxes, input boxes, and graphs
set(handles.windowType, 'Value', 1);
set(handles.blockSize, 'Value', 5);
set(handles.liveDataSampleRate, 'Value', 1);
set(handles.wavOutput, 'String', '');
set(handles.percentOverlap, 'Value', 1);
set(handles.avgpeakFreqA, 'String', '');
set(handles.avgpeakFreqB, 'String', '');
set(handles.avgpeakMagA, 'String', '');
set(handles.avgpeakMagB, 'String', '');
set(handles.instpeakFreqA, 'String', '');
set(handles.instpeakFreqB, 'String', '');
set(handles.instpeakMagA, 'String', '');
set(handles.instpeakMagB, 'String', '');
set(handles.expFilterResponse, 'Value', 1);
cla(handles.timeData);
cla(handles.instSpectrum);
cla(handles.avgSpectrum);

% --- Executes on button press in analysis.
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%initialize stop to off (0)
set(handles.stop, 'UserData', 0);
%grab data from buttons after being pushed to determine if live data or wav
%file read was pushed first
x = cell2mat(get(handles.analysis, 'Value'));
%decide which button was pushed and run based on parameters
%if wav analysis was pushed
if x(2) == 2
    %read in wav file
    [FileName, canc] = uigetfile('*.wav', 'Select the WAV file');
    %error if no wav file selected or canceled
    if canc == 0
        errordlg('No wav file was selected. ABORT!');
    end
    %decode file for wav file specifics
    [CalDate, SN, CalA, CalB, returnVal] = wavFileDecoder(FileName);
    %convert SN to string so that it can be displayed easier later if there
    %is no SN
    SN = num2str(SN);
    %determine if there is calibration information
    if returnVal == 0
        %ask for specific calibration information
        sensitivityPrompt = {sprintf('The sensor you are using does not contain calibration information. Please input the following in order to proceed with the analysis (nominal values shown).\n\nChannel A sensitivity (counts/(m/s^2)): '), 'Channel B sensitivity (counts/(m/s^2)): '};
        %set up the prompt specifics
        dlg_title = 'Sensitivity Input';
        num_lines = 1;
        %initialize the sensitivities
        def = {'33000', '65000'};
        sensitivity = inputdlg(sensitivityPrompt, dlg_title, num_lines, def);    %detect cancel
        %convert from string to double
        [CalA, CalB] = sensitivity{:};
        CalA = str2double(CalA);
        CalB = str2double(CalB);
        %get out of loop and continue callback
        returnVal = 1;
    end
    %determine if the wav file is properly formatted
    if  returnVal == -1
        errordlg('This wav file is not properly formatted');
    else
        %the wav file has Cal Data and is properly formatted
        [Y, sampleRate, nbits, ~] = wavread(FileName);
    end
%if live data was pushed    
elseif x(1) == 2
    %set up structure for recording
    continuousMode = true;
    hDSPAR = dsp.AudioRecorder('DeviceName', 'ASIO4ALL v2');
    hDSPAR.NumChannels = 2;
    hDSPAR.DeviceDataType = '24-bit integer';
    %decode the information on the sensor
    [SN, CalA, CalB, CalDate, version] = DigiDecoder();
    %determine if there is a cal date
    if version == 0
        CalDate = 'Not Available';
    end
        SN = num2str(SN);
        CalDate = char(CalDate);
end

%bring in user input specifics 
expFilterResponseArray = [0.25, 0.5, 0.75];
expFilterResponse = expFilterResponseArray(get(handles.expFilterResponse, 'Value'));
percentOverlapArray = [1, 0.25, 0.5, 0.75];
percentOverlap = percentOverlapArray(get(handles.percentOverlap, 'Value'));
blockSize = 2^(5+get(handles.blockSize, 'Value'));
windowTypeArray = {'Flattop', 'Hanning', 'Hann', 'Hamming', 'Blackman-harris'};
windowType = char(windowTypeArray(get(handles.windowType, 'Value')));
%decide specifics regarding live data or wav data
if x(1) == 2
    sampleRateArray = [8000, 11025, 16000, 22050, 32000, 44100, 48000];
    sampleRate = sampleRateArray(get(handles.liveDataSampleRate, 'Value'));
elseif x(2) == 2
    siz = wavread(FileName, 'size');
    secLength = siz(1)/sampleRate;
end
%initialize structure variables
hDSPAR.BufferSizeSource = 'Property';
hDSPAR.SampleRate = sampleRate;
hDSPAR.SamplesPerFrame = blockSize;
hDSPAR.BufferSize = 8192;
hDSPAR.QueueDuration = 5;
sample = 1;
%create axes for graphing
xtaxis = (1/hDSPAR.SampleRate)*(0:hDSPAR.SamplesPerFrame-1);
xfaxis = (hDSPAR.SampleRate/(hDSPAR.SamplesPerFrame))*(0:(hDSPAR.SamplesPerFrame/2)-1);
%initialize count
numofBlocks = 1;
gcf;
%continue loop until the stop button is pushed
while get(handles.stop,'UserData') ~= 1
    %% Frequency display
    if x(1) == 2
        %reads in the live data
        y = step(hDSPAR);
        if ~continuousMode
            release(hDSPAR);
        end
    elseif x(2) == 2
        %reads in wav data
        y = wavread(FileName, [sample, sample + blockSize-1]);
        %Select frequency display
    end
    %uses the sensitivity to convert to readable units
    yA = (2^23)/(9.90665*CalA)*y(:,1);
    yB = (2^23)/(9.90665*CalB)*y(:,2);
    %sets to timeData axes
    axes(handles.timeData);   
    %creates wav plot 
    if x(2) == 2
    p1(1) = plot(yA(1:blockSize),'r');
    p1(2) = plot(yB(1:blockSize), 'k');
    xlabel('Sample #');    
    xlim([0 blockSize]);
    ylim([-1 1]);
    grid on
    end
    %first time only
    if numofBlocks == 1
        %creates live plot for first time
        if x(1) == 2
            p1 = plot(xtaxis, yA,'r', xtaxis, yB, 'k');
            xlabel('Time (s)');
            xlim([0 max(xtaxis)]);
            grid on
        end
    else
        if x(1) == 2
            set(p1(1), 'XData', xtaxis, 'YData', yA);
            set(p1(2), 'Xdata', xtaxis, 'YData', yB);
            grid on
        end
    end
    title('Wave Function');
    ylabel(sprintf('Acceleration\n(g''s)'));
    %get time data for channel A
    % Compute spectrum for channel A
    xs1=spectralcalc(yA,1,hDSPAR.SamplesPerFrame-1,windowType); % scaling the halved channel to get the correct vibration amplitude
    % get time data for channel B
    % Compute spectrum for channel B
    xs2=spectralcalc(yB,1,hDSPAR.SamplesPerFrame-1,windowType);
    % averaging
    %sets to new graph
    axes(handles.instSpectrum);
    %determines dB
    xs1Mag = xs1.Magnitude(:)';
    xs2Mag = xs2.Magnitude(:)';
    if numofBlocks == 1
        % Channel average for first sample is itself
        xavgsum1 = xs1Mag;
        xavgsum2 = xs2Mag;
        %math and graph
        g2 = semilogx(xfaxis,20.*log10(xs1Mag),xfaxis,20.*log10(xs2Mag));
        xlabel('Frequency (Hz)');
        ylabel(sprintf('Instantaneous Acceleration\n(dB of fullscale)\nAvg #: %i',numofBlocks));
        grid on;
        xlim([0 sampleRate/2]);
        ylim([-160 60]);
        
    else
        %graphs after first time
        set(g2(2),'XData',xfaxis,'YData',20.*log10(xs1.Magnitude));
        set(g2(1),'XData',xfaxis,'YData',20.*log10(xs2.Magnitude));
        ylabel(handles.instSpectrum, sprintf('Instantaneous Acceleration\n(dB g)\nBlock #: %i',numofBlocks));
        title(handles.instSpectrum, 'Instantaneous Frequency Spectrum');
        grid on
    end
    %Compute peak statistics to adjust display
    %finds the top dB peak
    [~,peakFreqInd] = max(xs2(1).Magnitude);
    peakFreq = xfaxis(peakFreqInd);
    peakFreqMag = xs2(1).Magnitude(peakFreqInd);
    %sets to new graph, avgSpectrum
    axes(handles.avgSpectrum);
    %averages
    xavgsum1 = (expFilterResponse*xs1Mag+(1-expFilterResponse)*xavgsum1);
    xavgsum2 = (expFilterResponse*xs2Mag+(1-expFilterResponse)*xavgsum2);
    %first time avg Spectrum
    if numofBlocks == 1
        %determines average based on exponential filter response
        p2 = semilogx(xfaxis,20.*log10(xavgsum1),xfaxis,20.*log10(xavgsum2));
        xlabel('Frequency (Hz)');
        ylabel(sprintf('Avg Acceleration\n(dB of fullscale)\nAvg #: %i',numofBlocks));
        title(handles.avgSpectrum, 'Average Frequency Spectrum');
        grid on;
        xlim([0 sampleRate/2]);
        ylim([-160 60]);
    else
        %graphs after first time
        set(p2(2),'XData',xfaxis,'YData',20.*log10(xavgsum2));
        set(p2(1),'XData',xfaxis,'YData',20.*log10(xavgsum1));
        ylabel(handles.avgSpectrum, sprintf('Average Acceleration\n(dB g)\nAvg #: %i',numofBlocks));
        grid on
    end
    %in order to update figure and execute callbacks
    drawnow();
    %determine overlap for wav data
    if percentOverlap == 1
        %if overlap is 0, make sure there is no overlap
        jump = blockSize*percentOverlap;
        sample = sample+jump;
    else
        jump = blockSize*(1-percentOverlap);
        sample = sample + jump;
    end
     
    %determine max freq and dB for average
    [~,avgPeakFreqInd1] = max(xavgsum1);
    [~,avgPeakFreqInd2] = max(xavgsum2);
    avgpeakFreqA = xfaxis(avgPeakFreqInd1);
    avgpeakFreqB = xfaxis(avgPeakFreqInd2);
    avgpeakMagA = 20*log10(max(xavgsum1));
    avgpeakMagB = 20*log10(max(xavgsum2));
    
    %displays in text boxes in GUI
    set(handles.avgpeakFreqA, 'String', avgpeakFreqA);
    set(handles.avgpeakFreqB, 'String', avgpeakFreqB);
    set(handles.avgpeakMagA, 'String', avgpeakMagA);
    set(handles.avgpeakMagB, 'String', avgpeakMagB);
    
    %dtermines max freq and dB for instantaneous 
    [~,instPeakFreqInd1] = max(xs1.Magnitude);
    [~,instPeakFreqInd2] = max(xs2.Magnitude);
    instpeakFreqA = xfaxis(instPeakFreqInd1);
    instpeakFreqB = xfaxis(instPeakFreqInd2);
    instpeakMagA = 20*log10(max(xs1.Magnitude));
    instpeakMagB = 20*log10(max(xs2.Magnitude));
    
    %displays in text boxes in GUI
    set(handles.instpeakFreqA, 'String', instpeakFreqA);
    set(handles.instpeakFreqB, 'String', instpeakFreqB);
    set(handles.instpeakMagA, 'String', instpeakMagA);
    set(handles.instpeakMagB, 'String', instpeakMagB);
    
    %print out results accordingly
    if x(1) == 2
        if version == 1
            message = sprintf('Serial Number: %s\nCalibration Date: %s\nSample Rate (Hz): %i\nChannel A Sensitivity (counts/(m/s^2)): %i\nChannel B Sensitivity (counts/(m/s^2)): %i\nNumber of Blocks: %i\n', SN, CalDate, sampleRate, CalA, CalB, numofBlocks-1);
        end
    elseif x(2) == 2
        message = sprintf('WAV File Read: %s\nSerial Number: %s\nCalibration Date: %s\nSample Rate (Hz): %i\nNumber of Bits: %i\nLength of Time (s): %f\nChannel A Sensitivity (counts/(m/s^2)): %i\nChannel B Sensitivity (counts/(m/s^2)): %i\nNumber of Blocks: \n', FileName, SN, CalDate, sampleRate, nbits, secLength, CalA, CalB);
    end
    set(handles.wavOutput, 'String', message);
    %determine if there is more data from wav file to fit block size
    if x(2) == 2
        %if no more data, exit
        if (sample + blockSize - 1) >= siz(1)
            message = sprintf('WAV File Read: %s\nSerial Number: %s\nCalibration Date: %s\nSample Rate (Hz): %i\nNumber of Bits: %i\nLength of Time (s): %f\nChannel A Sensitivity (counts/(m/s^2)): %i\nChannel B Sensitivity (counts/(m/s^2)): %i\nNumber of Blocks: %i\n', FileName, SN, CalDate, sampleRate, nbits, secLength, CalA, CalB, numofBlocks-1);
            set(handles.wavOutput, 'String', message);
            return
        end
    end
    numofBlocks = numofBlocks +1;
    pause(0.000001)
end




function wavOutput_Callback(hObject, eventdata, handles)
% % hObject    handle to wavOutput (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of wavOutput as text
%        str2double(get(hObject,'String')) returns contents of wavOutput as a double
% Sample_Rate = get(hObject, 'String');

% --- Executes during object creation, after setting all properties.
function wavOutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function timeData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate timeData


% --- Executes on selection change in percentOverlap.
function percentOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to percentOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns percentOverlap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from percentOverlap

% --- Executes during object creation, after setting all properties.
function percentOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percentOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set options for percent overlap
set(hObject, 'String', {'0', '25', '50', '75'});



function sensitivityA_Callback(hObject, eventdata, handles)
% hObject    handle to sensitivityA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensitivityA as text
%        str2double(get(hObject,'String')) returns contents of sensitivityA as a double


% --- Executes during object creation, after setting all properties.
function sensitivityA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensitivityA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function instpeakMagA_Callback(hObject, eventdata, handles)
% hObject    handle to instpeakMagA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of instpeakMagA as text
%        str2double(get(hObject,'String')) returns contents of instpeakMagA as a double


% --- Executes during object creation, after setting all properties.
function instpeakMagA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instpeakMagA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgpeakFreqA_Callback(hObject, eventdata, handles)
% hObject    handle to avgpeakFreqA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgpeakFreqA as text
%        str2double(get(hObject,'String')) returns contents of avgpeakFreqA as a double


% --- Executes during object creation, after setting all properties.
function avgpeakFreqA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgpeakFreqA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function instpeakMagB_Callback(hObject, eventdata, handles)
% hObject    handle to instpeakMagB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of instpeakMagB as text
%        str2double(get(hObject,'String')) returns contents of instpeakMagB as a double


% --- Executes during object creation, after setting all properties.
function instpeakMagB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instpeakMagB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function instpeakFreqB_Callback(hObject, eventdata, handles)
% hObject    handle to instpeakFreqB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of instpeakFreqB as text
%        str2double(get(hObject,'String')) returns contents of instpeakFreqB as a double


% --- Executes during object creation, after setting all properties.
function instpeakFreqB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instpeakFreqB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function instpeakFreqA_Callback(hObject, eventdata, handles)
% hObject    handle to instpeakFreqA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of instpeakFreqA as text
%        str2double(get(hObject,'String')) returns contents of instpeakFreqA as a double


% --- Executes during object creation, after setting all properties.
function instpeakFreqA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instpeakFreqA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgpeakMagA_Callback(hObject, eventdata, handles)
% hObject    handle to avgpeakMagA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgpeakMagA as text
%        str2double(get(hObject,'String')) returns contents of avgpeakMagA as a double


% --- Executes during object creation, after setting all properties.
function avgpeakMagA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgpeakMagA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensitivityB_Callback(hObject, eventdata, handles)
% hObject    handle to sensitivityB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensitivityB as text
%        str2double(get(hObject,'String')) returns contents of sensitivityB as a double


% --- Executes during object creation, after setting all properties.
function sensitivityB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensitivityB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgpeakFreqB_Callback(hObject, eventdata, handles)
% hObject    handle to avgpeakFreqB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgpeakFreqB as text
%        str2double(get(hObject,'String')) returns contents of avgpeakFreqB as a double


% --- Executes during object creation, after setting all properties.
function avgpeakFreqB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgpeakFreqB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avgpeakMagB_Callback(hObject, eventdata, handles)
% hObject    handle to avgpeakMagB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avgpeakMagB as text
%        str2double(get(hObject,'String')) returns contents of avgpeakMagB as a double


% --- Executes during object creation, after setting all properties.
function avgpeakMagB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avgpeakMagB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in liveDataSampleRate.
function liveDataSampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to liveDataSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns liveDataSampleRate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from liveDataSampleRate


% --- Executes during object creation, after setting all properties.
function liveDataSampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to liveDataSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set up options for live data sample rate
set(hObject, 'String', {'8000', '11025', '16000', '22050', '32000', '44100', '48000'});
set(hObject, 'Value', 1);


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if clicked, set stop equal to 1 and stop the analyze callback
set(handles.stop, 'UserData', 1);
guidata(hObject, handles);


% --- Executes on selection change in expFilterResponse.
function expFilterResponse_Callback(hObject, eventdata, handles)
% hObject    handle to expFilterResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns expFilterResponse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from expFilterResponse


% --- Executes during object creation, after setting all properties.
function expFilterResponse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expFilterResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set up options for the exponential filter response and initialize to slow
set(hObject, 'String', {'Slow', 'Medium', 'Fast'});
set(hObject, 'Value', 1);
