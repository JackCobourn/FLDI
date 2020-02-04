function NIDAQ_Scope(DAQSession)
% NIDAQ_Scope Use the NI DAQ as an osilloscope
plotter = @(src,event) plot(event.TimeStamps, event.Data);

lh = addlistener(DAQSession,'DataAvailable', plotter); 
DAQSession.IsContinuous = true;

startBackground(DAQSession);



function plotDAQ(src,event)
     plot(event.TimeStamps, event.Data)
end
end

%DAQSession.Stop to stop