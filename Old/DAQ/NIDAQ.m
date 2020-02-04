% NI DAQ Session Connector
function [DAQSession,DAQ,AIn,FLDIChannel] = NIDAQ(~)
    
    %Function to connnect the NI DAQ currently connected to the PC,
    %configure Analog Input 0 to +1 to -1 Volts
    %sample rate to the max, 1 sec aquisiition.
    %and the return all needed structures.




    daqreset; %Reset Matlab DAQ controller
    
    %Get info on connected devices, if theres more than one youre gonna have to specify that
    DAQ = daq.getDevices; 
    %get the devices subsystems
    subsys = DAQ.Subsystems;
    
    %get the Analog Input subsys
    AIn = subsys(1);
    
    %create a DAQ session for all NI DAQs
    DAQSession = daq.createSession('ni'); 
    
    %tell matlab about the channel you are using using the device ID and
    %daq session. AIn.ChannelNames{1,1} is AI0. Also need to tell it youre measuring Voltage.
    FLDIChannel = addAnalogInputChannel(DAQSession,DAQ.ID, AIn.ChannelNames{1,1}, 'Voltage');
    
    %set the aquisition rate to max bc ¯\_(?)_/¯
    DAQSession.Rate = DAQSession.RateLimit(2);
    
    %set the range to +1 to -1, which is the minimum
    %(1,2) is +2,-2. (1,3) is +5,-5, (1,4) is +10,-10.
    FLDIChannel.Range = subsys(1).RangesAvailable(1,1);
    
end