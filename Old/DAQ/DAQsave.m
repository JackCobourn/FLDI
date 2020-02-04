function [Data,Times] = DAQsave(src,event)
    Data = event.Data;
    Times = event.TimeStamps;
end