%% Oscope v2
%function VisaInterface = OscopeVisa(~)

%connect to the oscope with all of: matlab connection to NI Visa, NI Visa,
%and Lecroy connection to Ni Visa

%Previous to using this script you need to create a logical name for your
%hardware and driver compbo in NI MAX

%then use makemid() to create a matlab wraper for the lacroy driver that
%has been interpredted by NI MAX

%this allows for both high level object calls though visual basic (ie"vbs
%'app.' " calls to the machine
%and low level "ALL_CAPS" calls to the old IEEE488.2 Commands (they don't
%actually need to be all caps)
%see maui-remote-control-and-automation-manual.pdf from lecroy Wavesurfer
%510 resources

%however, the driver I edited, to create a matlab device object, only makes
%the low level calls.

VisaInterface = visa('ni','TCPIP0::192.168.137.2::INSTR')
deviceObj = icdevice('lecroy_WaveSurvfer510_VISA.mdd', interfaceObj);

fopen(VisaInterface)
VisaInterface.timeout=5

%use a low level call to turn off header and a high level to set the whole
%device to defult
fwrite(VisaInterface,"COMM_HEADER OFF");
fwrite(VisaInterface,"vbs 'app.settodefaultsetup'");

%combine channels 1-2 and 3-4 in order to increase availabe mem per
%channel, as we are only using channel 3 for data and channel 1 for trigger
%(maybe in M4 and M2, not on the benchtop)
fprintf(VisaInterface, "COMBINE_CHANNELS 2");
%set mem to max
fprintf(VisaInterface, "MEMORY_SIZE 20E+6");

%Calculate the largest needed buffer for waveform transfer
MAX_HEADER = 1e3;
memorySize = str2num(query(VisaInterface, 'Memory_SIZe?'));
bufferSize = 2*(memorySize+MAX_HEADER);


%close the connection, change the buffers, and re open
fclose(VisaInterface)
VisaInterface.InputBufferSize = bufferSize;
VisaInterface.OutputBufferSize = bufferSize;
fopen(VisaInterface)

sampleElement = cast(1, 'int8');
sampleElementProperties = whos('sampleElement');
elementSize = sampleElementProperties.bytes;
elementsToRead = get(VisaInterface, 'InputBufferSize')/elementSize; 

query(VisaInterface, 'CFMT?')

%fwrite(VisaInterface,"VBS 'app.Memory.M1.Out.Result.DataArray'")
fwrite(VisaInterface,"VBS? 'return=app.acquisition.C3.Out.Result.DataArray'")
fwrite(VisaInterface,"VBS? 'return=app.acquisition.C3.Out.Result.HorizontalOffset'")

%fwrite(VisaInterface, "M1:WF? ALL")
%fwrite(VisaInterface, "M1:WF? DESC")
fwrite(VisaInterface, "M1:WF? DAT1")
msg = [];
data = [];
while isempty(msg) %(~strcmp(msg, foundEoiMsg) || isempty(findstr(msg, timeoutSubMsg)))
        [newData, count, msg] = fread(VisaInterface, elementsToRead, 'int8');
        data = [data; newData];
end


%Check out what the waveform template looks like
%LeCroy scopes publish the format of the waveform (WAVEDESC) with the 
% TEMPLATE? query.  By reading and parsing that, the code can determine
% where all of the fields in the waveform are regardless of changes to
% the scope.  
template = query(VisaInterface, 'TEMPLATE?');

%% Test Memory Time
tic
[y,~,info]= invoke(WVE,'readwaveform','C3',false);
toc
tic
invoke(WVE,'storewaveform','C3','M2');
[y2,~,info2] = invoke(WVE,'readwaveform','M2',false);
toc
%% Learn to proceess template
idxWaveDesc = regexp(template,'(WAVEDESC: BLOCK)|(ENDBLOCK)','start'); %find the index where the actual description starts

templateStruct=regexp(template(idxWaveDesc(1): idxWaveDesc(2)),... %find the offsets, info, and data types
    '<(?<idxStart>[0-9 ]*)>[\s]*(?<name>[\w]+):\s*(?<type>[\w]+)','names');

% Convert the indexes (currently character arrays) to numbers and account
% for the initial byte offset (headerOffset) to index into the waveform.  
% Additionally, convert the known types in the waveform to MATLAB types.  
% Not all possible types are accounted for.

headerOffset = 16; 

for i=1:length(templateStruct)
    % convert to numbers (from chars) and add byte offset
    templateStruct(i).idxStart = ...
        str2num(templateStruct(i).idxStart) + headerOffset;
    
    % convert template type names to MATLAB types
    templateValues = {'string', 'byte', 'word',  'long',  'float',  'double', 'text', 'unit_definition'};
    matlabValues = {  'char',   'int8', 'int16', 'int32', 'single', 'double', 'uint8', 'uint8'};
    idx = strmatch(templateStruct(i).type, templateValues);
    if ~isempty(idx)
        templateStruct(i).type = matlabValues{idx};
    end
end

% This parses all but the last item of the of the template because
% the last item isn't necessary and the code is cleaner.
% Store the starting and ending field indices under fieldnames (where the 
% name is the property name from the waveform header) of a structure.
for i=1:length(templateStruct)-1
    wd.(templateStruct(i).name).idxStart = ...
        templateStruct(i).idxStart;
    wd.(templateStruct(i).name).idxEnd = ...
        templateStruct(i+1).idxStart - 1; %it ends at the start of the next item
    wd.(templateStruct(i).name).type = templateStruct(i).type;
end

%%%%%%%%%%%%%%%%%%%%%%
%%TEST DRIVER
% initializes everything to empty/false/unknown
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
 
%checks to see if the source input is a valid value
validValues = {'C1', 'C2', 'C3', 'C4', 'TA', 'TB', 'TC', 'TD', 'M1', 'M2', 'M3', 'M4'};
scopeValues = {'C1', 'C2', 'C3', 'C4', 'TA', 'TB', 'TC', 'TD', 'M1', 'M2', 'M3', 'M4'};
idx = strmatch(source,validValues, 'exact');

%gives error if the channel is not a valid value
if (isempty(idx))
    error('Invalid SOURCE. CHANNEL must be one of: C1, C2, C3, C4, TA, TB, TC, TD, M1, M2, M3, M4');
end

%sets the source to its property value
trueSource = scopeValues{idx};

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
        %set(obj, 'Precision', oldPrecision);
        %set(obj, 'ByteOrder', oldByteOrder);
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
         header = raw(2:end)'; %deals with the fact that All is shorter than DESC, which is what the template parser expects, andwhich I don't want to rewrite.
         
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
         %set(obj, 'ByteOrder', oldByteOrder);
         error('An error occurred while reading the waveform.');
     end
else
  error('Ascii precsision not allowed for WF transfer');
end


xunit = wd.getHORUNIT(header);
yunit = wd.getVERTUNIT(header);
ygain = wd.getVERTICAL_GAIN(header);
yintercept = wd.getVERTICAL_OFFSET(header);
length = wd.getWAVE_ARRAY_COUNT(header);
xinterval = wd.getHORIZ_INTERVAL(header);
xintercept = wd.getHORIZ_OFFSET(header);
Info = struct('xunit',xunit,'yunit',yunit,'ygain',ygain,'yintercept',yintercept,'length',length,'xinterval',xinterval,'xintercept',xintercept);
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
%set(obj, 'Precision', oldPrecision);
%set(obj, 'ByteOrder', oldByteOrder);r);