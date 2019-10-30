function [VisaInterface, Wavesurfer510, Waveform] = OscopeVisaConnect_V2(~)
instrreset
clear VisaInterface Wavesurfer510 Waveform Wave
VisaInterface = visa('ni','TCPIP0::192.168.137.2::INSTR');
Wavesurfer510 = icdevice('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf\lecroy_WaveSurvfer510_VISA.mdd', VisaInterface);
%Wavesurfer510 = icdevice('lecroy_8600a.mdd', VisaInterface);

connect(Wavesurfer510)
Waveform=Wavesurfer510.Waveform;

VisaInterface
Wavesurfer510
Waveform
end


%%%%%%%%%%%%%%%%%%%%%TEST diver 
%Set stuff by hand on Oscope, its easier for now
%set(Waveform,'Precision','int8')
%[Y1,X1,YUNIT1,XUNIT1,HEADER1] = invoke(Waveform,'readwaveform','channel3');
%set(Waveform,'Precision','int16');
y= invoke(Waveform,'readwaveform','C3');
%set(Waveform,'Precision','int16')
%[Y2,X2,YUNIT2,XUNIT2,HEADER2] = invoke(Waveform,'readwaveform','channel3');

%get(Wavesurfer510.Acquisition,'T
%template2 = query(VisaInterface, 'TEMPLATE?');

%%%%%%%%%%%%%%%%%%%%%TEST diver 

scale = true;
y = [];
x = [];
yunit = 'unknown';
xunit = 'unknown';
ygain = 'unknown';
yintercept = 'unknown';
length = 'unknown';
xinterval = 'unknown';
xintercept = 'unknown';
Info = struct();
 
% checks to see if the source input is a valid value
%validValues = {'C1', 'C2', 'C3', 'C4', 'TA', 'TB', 'TC', 'TD', 'M1', 'M2', 'M3', 'M4'};
%scopeValues = {'C1', 'C2', 'C3', 'C4', 'TA', 'TB', 'TC', 'TD', 'M1', 'M2', 'M3', 'M4'};
%idx = strmatch(source,validValues, 'exact');

% gives error if the channel is not a valid value
% if (isempty(idx))
%     error('Invalid SOURCE. CHANNEL must be one of: C1, C2, C3, C4, TA, TB, TC, TD, M1, M2, M3, M4');
% end

% sets the source to its property value
trueSource = source;%scopeValues{idx};

% checks to see if the scale entered is a valid value
if exist('scaleWaveform', 'var')
        scale = scaleWaveform;
end

device = get(obj, 'Parent');
util = get(device, 'Util');

% get the driverdata
driverData = get(device, 'DriverData');
wd = driverData.wd;

% stores the old user setings
oldPrecision = get(obj, 'Precision');
oldByteOrder = get(obj, 'ByteOrder');

% changes the precision and byteorder so that wave form can be read
%set(obj, 'Precision', 'int8');
set(obj, 'ByteOrder', 'littleEndian');

% gets wave from the instrument
if strcmp(oldPrecision,'int8')
    try 
        % Issue the curve transfer command.
        to_interface = sprintf('%s:WaveForm? ALL',trueSource);
        invoke(util, 'sendcommand', to_interface);

        % get wave from the instrument
        raw = invoke(util, 'readbin', 'int8');

        % separates the data from the descriptor block 
       
        header = raw(1:361)';
        y = raw(362:end)';
    catch
        % return the settings back to what the user had them and display
        % an error that tells why the waveform was not read
        %set(obj, 'Precision', oldPrecision);
        set(obj, 'ByteOrder', oldByteOrder);
        error(lasterr)  
    end

    % make sure that data was
    if (isempty(raw))
        % returns the settings back to what the user had them and displays
        % an error that the waveform did not get read in
        set(obj, 'Precision', oldPrecision);
        set(obj, 'ByteOrder', oldByteOrder);
        error('An error occurred while reading the waveform.');
    end
  elseif strcmp(oldPrecision,'int16')
     try 
         % Issue the Header transfer command.
         to_interface = sprintf('%s:WaveForm? DESC',trueSource);
         invoke(util, 'sendcommand', to_interface);
 
         % get wave from the instrument
         raw = invoke(util, 'readbin', 'int8');
 
         % separates the data from the descriptor block 
         header = raw(2:end)'; %deals with the fact that 
         
         to_interface = sprintf('%s:WaveForm? DAT1',trueSource);
         invoke(util, 'sendcommand', to_interface);
 
         % get wave from the instrument
         raw = invoke(util, 'readbin', 'int16');
                
         % separates the data from the descriptor block 
         y = raw(9:end)';
%          %then get the data
%          to_interface = sprintf('%s:WaveForm? WAVE_ARRAY_1',trueSource);
%          invoke(util, 'sendcommand', to_interface);
%  
%          % get wave from the instrument
%          raw = invoke(util, 'readbin', 'int16');        
%          
%          y = raw(8:end); %need to check this
 
         
     catch
         % return the settings back to what the user had them and display
         % an error that tells why the waveform was not read
         %set(obj, 'Precision', oldPrecision);
         set(obj, 'ByteOrder', oldByteOrder);
         error(lasterr)  
     end
 
     % make sure that data was
     if (isempty(raw))
         % returns the settings back to what the user had them and displays
         % an error that the waveform did not get read in
         %set(obj, 'Precision', oldPrecision);
         set(obj, 'ByteOrder', oldByteOrder);
         error('An error occurred while reading the waveform.');
     end
else
  error('Ascii precsision not allowed for WF transfer');
end


xunit = wd.getHORUNIT(header);
yunit = wd.getVERTUNIT(header);
ygain = wd.getVERTICAL_GAIN(header);
yintercept =wd.getVERTICAL_OFFSET(header);
length = wd.getWAVE_ARRAY_COUNT(header);
xinterval = wd.getHORIZ_INTERVAL(header);
xintercept = wd.getHORIZ_OFFSET(header);
Info = struct('xunit',xunit,'yunit',yunit,'ygain',ygain,'yintercept',yintercept,'length',length,'xinterval',xinterval,'xintercept',xintercept)
% Rescale the waveform is requested.
if scale
    % rescale the Y
    y = y .* wd.getVERTICAL_GAIN(header) - wd.getVERTICAL_OFFSET(header);
end

% Calculate the horizontal points.
x = ((0:double(wd.getWAVE_ARRAY_COUNT(header))).*wd.getHORIZ_INTERVAL(header)...
    + wd.getHORIZ_OFFSET(header))';

x = x';
y = y';

% Restore initial settings
set(obj, 'Precision', oldPrecision);
set(obj, 'ByteOrder', oldByteOrder);