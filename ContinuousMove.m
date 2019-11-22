%ContinuousMotion

% tic
% MotorY.moveto(0)
% toc
% data = [];
% saver = @(src,event,data) [data;event.Data];
prepare(DAQSession)


DAQSession.DurationInSeconds = 2;
DAQSession.Rate = 100e3;
lh = addlistener(DAQSession,'DataAvailable', @plotDAQ); 
startBackground(DAQSession);
MotorY.moveto(50);


DAQSession.IsNotifyWhenDataAvailableExceedsAuto =false;
DAQSession.NotifyWhenDataAvailableExceeds = DAQSession.NumberOfScans;

