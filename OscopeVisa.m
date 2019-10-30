%% Oscope v2
function VisaInterface = OscopeVisa(~)

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
%Check out what the waveform template looks like
%LeCroy scopes publish the format of the waveform (WAVEDESC) with the 
% TEMPLATE? query.  By reading and parsing that, the code can determine
% where all of the fields in the waveform are regardless of changes to
% the scope.  
template = query(VisaInterface, 'TEMPLATE?');

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