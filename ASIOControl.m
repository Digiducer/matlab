% Load ASIO Control
% Requires ASIOControl.ocx registered using regsvr32
% Tim Collins, Dec 2008.

function h = ASIOControl(varargin)
% Check for an already open control
figs = get(0,'children');
for n = figs'
    if(isa(get(n, 'userdata'),'COM.ASIOCONTROL_ASIOControl_1'))
        h = get(n, 'userdata');
        if(~isempty(varargin))
            if(strcmp(varargin{1},'Resize'))
                pos = get(n, 'position');
                h.move([0 0 pos(3:4)]);
            else
                openASIODevice(h, varargin{1})
            end    
        end
        return
    end
end
% If none open, start a new one
fig = figure;
p = get(fig,'position');
set(fig,'position',[p(1:2) 400 200]);
h = actxcontrol('ASIOCONTROL.ASIOControl.1',[0 0 400 200],fig);
set(fig,'ResizeFcn','ASIOControl Resize;','Toolbar','none','NextPlot','new','userdata',h);
mh = uimenu(fig, 'Label', 'ASIO');
uimenu(mh, 'Label','Properties', 'callback', 'propedit(get(gcf,''UserData''))');

% If arg(1) is present, this is the device name
if(~isempty(varargin))
    openASIODevice(h, varargin{1})
end

function openASIODevice(h, name)
h.ASIODeviceName = name;
% Default set-up is stereo-in, stereo out
% Properties can be changed later if needed
h.channelsInputting = [1 1];
h.channelsOutputting = [1 1];
